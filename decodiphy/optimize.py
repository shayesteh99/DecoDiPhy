import numpy as np
import scipy.sparse as sp
import scipy
import time

import osqp
import tracemalloc

class AnchorOSQP:
    def __init__(
        self,
        d, l_T, C_T, D_T, init_anchors,
        eps_abs=1e-8, eps_rel=1e-8, max_iter=200000,
        polish=True,
        # rho / adaptive rho
        rho=0.1, adaptive_rho=True,
        adaptive_rho_interval=25,
        adaptive_rho_tolerance=5.0,
        adaptive_rho_fraction=0.4,
        # scaling options
        scaling=False, scaled_termination=False,
        # regularizer params for numerical PSD
        reg_rel=1e-9, reg_abs=1e-12,
        verbose=False
    ):
        self.d = np.asarray(d, dtype=np.float64).reshape(-1)
        self.l_T = np.asarray(l_T, dtype=np.float64).reshape(-1)
        self.C_T = np.asarray(C_T, dtype=np.float64)
        self.D_T = np.asarray(D_T, dtype=np.float64)

        self.anchors = np.asarray(init_anchors, dtype=int).copy()
        self.k = len(self.anchors)
        self.L = len(self.d)
        self.n_vars = 2 * self.k + 1
        self.last_z = np.zeros(self.n_vars)

        # solver params
        self.eps_abs = float(eps_abs)
        self.eps_rel = float(eps_rel)
        self.max_iter = int(max_iter)
        self.polish = bool(polish)
        self.verbose = bool(verbose)

        # rho / adaptive rho
        self.rho = float(rho)
        self.adaptive_rho = bool(adaptive_rho)
        self.adaptive_rho_interval = int(adaptive_rho_interval)
        self.adaptive_rho_tolerance = float(adaptive_rho_tolerance)
        self.adaptive_rho_fraction = float(adaptive_rho_fraction)

        # scaling
        self.scaling = bool(scaling)
        self.scaled_termination = bool(scaled_termination)

        # regularization
        self.reg_rel = float(reg_rel)
        self.reg_abs = float(reg_abs)

        # OSQP state
        self.prob = None
        self.P_pattern = None
        self.initialized = False

        # ---------------- constraints ----------------
        row_data, col_data, val_data = [], [], []
        l_list, u_list = [], []
        r = 0

        # p >= 0
        for i in range(self.k):
            row_data.append(r)
            col_data.append(i)
            val_data.append(1.0)
            l_list.append(0.0)
            u_list.append(np.inf)
            r += 1

        # x >= 0
        for i in range(self.k):
            row_data.append(r)
            col_data.append(self.k + i)
            val_data.append(1.0)
            l_list.append(0.0)
            u_list.append(np.inf)
            r += 1

        # x <= p
        for i in range(self.k):
            row_data.append(r)
            col_data.append(i)
            val_data.append(-1.0)

            row_data.append(r)
            col_data.append(self.k + i)
            val_data.append(1.0)

            l_list.append(-np.inf)
            u_list.append(0.0)
            r += 1

        # y >= 0
        row_data.append(r)
        col_data.append(2 * self.k)
        val_data.append(1.0)
        l_list.append(0.0)
        u_list.append(np.inf)
        r += 1

        # sum(p) = 1
        for i in range(self.k):
            row_data.append(r)
            col_data.append(i)
            val_data.append(1.0)
        l_list.append(1.0)
        u_list.append(1.0)
        r += 1

        self.A = sp.csc_matrix(
            (val_data, (row_data, col_data)),
            shape=(r, self.n_vars)
        )
        self.l_vec = np.array(l_list, dtype=float)
        self.u_vec = np.array(u_list, dtype=float)

        # ---------------- M matrix ----------------
        self.M_cols = np.zeros((self.L, self.n_vars), dtype=np.float64)
        for t, a in enumerate(self.anchors):
            self.M_cols[:, t] = self.D_T[a].T
            self.M_cols[:, self.k + t] = self.C_T[a].T * self.l_T[a]
        self.M_cols[:, 2 * self.k] = 1.0

        M = self.M_cols

        self.G = self.M_cols.T @ self.M_cols 

        # ---------------- P initialization ----------------
        P_full = 2.0 * (M.T @ M)
        P_full = 0.5 * (P_full + P_full.T)

        max_diag = np.max(np.abs(np.diag(P_full))) if P_full.size > 0 else 0.0
        reg = max(self.reg_abs, self.reg_rel * max(1.0, max_diag))
        P_full[np.diag_indices_from(P_full)] += reg
        n = self.n_vars

        P_template = sp.triu(sp.csc_matrix(P_full), format="csc")

        indptr = P_template.indptr
        indices = P_template.indices

        idx_map = -np.ones((n, n), dtype=int)
        rows, cols = [], []
        pos = 0

        for col in range(n):
            for p in range(indptr[col], indptr[col + 1]):
                row = indices[p]
                idx_map[row, col] = pos
                idx_map[col, row] = pos
                rows.append(row)
                cols.append(col)
                pos += 1

        self.P_template = P_template
        self.P_data = P_template.data.copy()
        self.idx_map = idx_map

        self.q = -2.0 * (M.T @ self.d)

        # if SOLVER == "osqp":
        self.prob = osqp.OSQP()

        self.prob.setup(
            P=sp.csc_matrix((self.P_data, (rows, cols)), shape=(n, n)),
            q=self.q,
            A=self.A,
            l=self.l_vec,
            u=self.u_vec,
            eps_abs=self.eps_abs,
            eps_rel=self.eps_rel,
            max_iter=self.max_iter,
            polish=self.polish,
            rho=self.rho,
            adaptive_rho=self.adaptive_rho,
            adaptive_rho_interval=self.adaptive_rho_interval,
            adaptive_rho_tolerance=self.adaptive_rho_tolerance,
            adaptive_rho_fraction=self.adaptive_rho_fraction,
            scaling=self.scaling,
            scaled_termination=self.scaled_termination,
            verbose=self.verbose
            )

        self.initialized = True
        self.last_anchors = self.anchors.copy()

    # --------------------------------------------------

    def solve(self, anchors_new):
        tracemalloc.start()
        t0 = time.time()
        anchors_new = np.asarray(anchors_new, dtype=int).copy()

        diffs = np.where(self.last_anchors != anchors_new)[0]

        # CASE 0: no change
        if diffs.size == 0:
            self.prob.warm_start(x=self.last_z)
            res = self.prob.solve()
            self.last_z = res.x.copy()
            t1 = time.time()

        # CASE 1: single anchor change
        elif diffs.size == 1:
            pos = int(diffs[0])
            new_anchor = int(anchors_new[pos])

            # --- update M columns ---
            new_col_p = self.D_T[new_anchor].T
            new_col_x = self.C_T[new_anchor].T * self.l_T[new_anchor]

            self.M_cols[:, pos] = new_col_p
            self.M_cols[:, self.k + pos] = new_col_x

            self.anchors[pos] = new_anchor
            self.last_anchors = self.anchors.copy()

            # --- update q ---
            self.q[pos] = -2.0 * np.dot(new_col_p, self.d)
            self.q[self.k + pos] = -2.0 * np.dot(new_col_x, self.d)
            self.q[-1] = -2.0 * self.d.sum()

            # --- incremental G update ---
            cols_to_update = [pos, self.k + pos]

            for i in cols_to_update:
                col_i = self.M_cols[:, i]
                # update row i of G
                self.G[i, :] = self.M_cols.T @ col_i
                self.G[:, i] = self.G[i, :]  # symmetry

            # constant column (2k) never changes → G[:, 2k] unchanged

            # --- write G into sparse P_data ---
            # P = 2*G + reg*I
            for i in cols_to_update:
                for j in range(self.n_vars):
                    r, c = min(i, j), max(i, j)
                    idx = self.idx_map[r, c]
                    if idx >= 0:
                        val = 2.0 * self.G[i, j]
                        if i == j:
                            val += self.reg_abs
                        self.P_data[idx] = val

            # --- update OSQP ---
            t1 = time.time()
            self.prob.update(Px=self.P_data, q=self.q)
            t2 = time.time()
            self.prob.warm_start(x=self.last_z)
            res = self.prob.solve()
            self.last_z = res.x.copy()


        # CASE 2: multiple changes
        else:
            print("case 2 hit!!!")
            for t, a in enumerate(anchors_new):
                self.M_cols[:, t] = self.D_T[a].T
                self.M_cols[:, self.k + t] = self.C_T[a].T * self.l_T[a]
            self.M_cols[:, 2 * self.k] = 1.0

            self.anchors = anchors_new.copy()
            self.last_anchors = anchors_new.copy()

            P_full = 2.0 * (self.M_cols.T @ self.M_cols)
            P_full = 0.5 * (P_full + P_full.T)
            P_full[np.diag_indices_from(P_full)] += self.reg_abs

            P_csc = sp.triu(sp.csc_matrix(P_full), format="csc")
            self.P_data = P_csc.data.copy()
            self.q = -2.0 * (self.M_cols.T @ self.d)

            self.prob.update(Px=self.P_data, q=self.q)
            self.prob.warm_start(x=self.last_z)
            res = self.prob.solve()
            self.last_z = res.x.copy()

        if res.info.status_val not in (1, 2):
            return float("inf"), None, None, None, t2 - t1, time.time() - t2

        z = res.x
        p = z[:self.k]
        x = z[self.k:2 * self.k]
        y = z[-1]

        DA = self.M_cols[:, :self.k]
        CAl = self.M_cols[:, self.k:2 * self.k]
        obj = np.sum((self.d - (DA @ p + CAl @ x + y)) ** 2)

        ratio = np.zeros(self.k)
        nz = p != 0
        ratio[nz] = np.minimum(1.0, np.abs(x[nz] / p[nz]))

        
        current, peak = tracemalloc.get_traced_memory()
        # print(f"Current memory: {current / 1e6:.2f} MB")
        # print(f"Peak memory: {peak / 1e6:.2f} MB")

        tracemalloc.stop()
        return obj, p, ratio, y, time.time() - t0, peak / 1e6

def get_optimal_obj(problem, d, l_T, C_T, D_T, anchors, k):
    t = time.time()
    # DA = D[:, anchors]
    DA = D_T[anchors].T
    # CAl = C[:, anchors] * l[anchors].T
    CAl = C_T[anchors].T * l_T[anchors]

    # n = len(index_to_node)
    L = len(d)

    # p = cp.Variable(k)   
    # x = cp.Variable(k)   
    # y = cp.Variable()

    problem['DA'].value = DA
    problem['CAl'].value = CAl
    problem['d'].value = d

    # term1 = DA @ p
    # term2 = CAl @ x
    # term3 = np.ones(L) * y

    # residuals = d - (term1 + term2 + term3)
    # objective = cp.Minimize(cp.sum_squares(residuals))

    # constraints = [
    # p >= 0,
    # x >= 0,
    # x <= p,
    # y >= 0,
    # cp.sum(p) == 1,
    # ]

    # problem = cp.Problem(objective, constraints)
    # print(problem.size_metrics)
    # print(problem.size_metrics.__dict__)
    try:
        problem['problem'].solve(verbose=False)
    except cp.error.SolverError:
        print("Solver reached max iterations. Using the last available solution.")
        e = time.time()
        return float('Inf'), None, None, None, e-t

    e = time.time()
    return problem['problem'].value, problem['p'].value, np.array([min(i, 1) for i in np.fabs(problem['x'].value/problem['p'].value)]), problem['y'].value, e-t, e-t

def build_cvxpy_problem(L, k):
    p = cp.Variable(k)
    x = cp.Variable(k)
    y = cp.Variable()

    DA_param = cp.Parameter((L, k))
    CAl_param = cp.Parameter((L, k))
    d_param = cp.Parameter(L)

    residuals = d_param - (DA_param @ p + CAl_param @ x + y)

    objective = cp.Minimize(cp.sum_squares(residuals))

    constraints = [
        p >= 0,
        x >= 0,
        x <= p,
        y >= 0,
        cp.sum(p) == 1,
    ]

    problem = cp.Problem(objective, constraints)

    res= {}
    res['problem'] = problem
    res['p'] = p
    res['x'] = x
    res['y'] = y
    res['DA'] = DA_param
    res['CAl'] = CAl_param
    res['d'] = d_param

    return res


class AnchorCVXPY:
    def __init__(self, d, l_T, C_T, D_T, init_anchors):
        self.k = len(init_anchors)
        self.d = d
        self.l_T = l_T
        self.C_T = C_T
        self.D_T = D_T

        self.problem = build_cvxpy_problem(len(d), self.k)


    def solve(self, anchors_new):
        obj, p, x, y, t = get_optimal_obj(self.problem, self.d, self.l_T, self.C_T, self.D_T, anchors_new, self.k)
        return obj, p, x, y, t, 0.0




def solve_osqp_once(d, l_T, C_T, D_T, anchors, index_to_node, k):
    t = time.time()
    DA = D_T[anchors].T
    CAl = C_T[anchors].T * l_T[anchors]
    L = len(d)

    M = np.hstack([DA, CAl, np.ones((L, 1))])
    P = 2.0 * (M.T @ M)
    q = -2.0 * (M.T @ d)

    n_vars = 2 * k + 1
    row_data, col_data, val_data = [], [], []
    l_list, u_list = [], []
    r = 0

    for i in range(k):
        row_data.append(r); col_data.append(i); val_data.append(1.0)
        l_list.append(0.0); u_list.append(np.inf)
        r += 1

    for i in range(k):
        row_data.append(r); col_data.append(k + i); val_data.append(1.0)
        l_list.append(0.0); u_list.append(np.inf)
        r += 1

    for i in range(k):
        row_data.append(r); col_data.append(i); val_data.append(-1.0)
        row_data.append(r); col_data.append(k + i); val_data.append(1.0)
        l_list.append(-np.inf); u_list.append(0.0)
        r += 1

    row_data.append(r); col_data.append(2 * k); val_data.append(1.0)
    l_list.append(0.0); u_list.append(np.inf)
    r += 1

    for i in range(k):
        row_data.append(r); col_data.append(i); val_data.append(1.0)
    l_list.append(1.0); u_list.append(1.0)
    r += 1

    A = sp.csc_matrix((val_data, (row_data, col_data)), shape=(r, n_vars))
    l_vec = np.array(l_list, dtype=float)
    u_vec = np.array(u_list, dtype=float)
    P_csc = sp.csc_matrix(P)

    prob = osqp.OSQP()
    prob.setup(
        P=P_csc, q=q, A=A, l=l_vec, u=u_vec,
        eps_abs=1e-8, eps_rel=1e-8, max_iter=200000,
        polish=True, rho=0.1, adaptive_rho=True,
        adaptive_rho_interval=25, adaptive_rho_tolerance=5.0,
        adaptive_rho_fraction=0.4, scaling=False,
        scaled_termination=False, verbose=False,
    )
    res = prob.solve()

    if res.info.status_val not in [1, 2]:
        print("OSQP did not converge")
        return float('inf'), None, None, None, time.time() - t

    z_opt = res.x
    p_opt = z_opt[:k]
    x_opt = z_opt[k:2 * k]
    y_opt = z_opt[-1]
    obj_val = np.sum((d - (DA @ p_opt + CAl @ x_opt + y_opt)) ** 2)

    return obj_val, p_opt, np.array([min(i, 1) for i in np.fabs(x_opt / p_opt)]), y_opt, time.time() - t
