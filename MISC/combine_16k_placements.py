import sys
import os
import argparse
import time
import numpy as np
from treeswift import *
import treeswift
import random
import json
# from kneed import KneeLocator

def choose_round(all_rounds, total_p, threshold):
	#knee 
	losses=[]
	ks = []
	for a in all_rounds:
		if 'p_value' in a:
			losses.append(np.log(a['loss']))
			ks.append(a['k'])

	if len(ks) == 1:
		knee= 1
	else:
		knee = KneeLocator(ks, losses, curve="convex", direction="decreasing", S = total_p * threshold).elbow
	# print(knee)

	for a in all_rounds:
		if 'p_value' in a and a['k'] == knee:
			return a
	return all_rounds[-1]


	# max_k * threshold
	max_k = 0
	for a in all_rounds:
		if 'p_value' in a:
			max_k = a['k']

	for a in all_rounds:
		if 'p_value' in a:
			if a['k'] >= max_k * threshold:
				return a
	return all_rounds[-1]


	# loss diff * p
	prev_loss = None
	for a in all_rounds:
		if 'p_value' in a:
			if not prev_loss:
				prev_loss = a['loss']
				continue
			diff = total_p * (prev_loss - a['loss']) / prev_loss
			if diff < threshold:
				return a


			# if a['loss'] * total_p < threshold:
	# print(ks, losses)


	return all_rounds[-1]

def choose_round_diff(all_rounds, total_p, threshold):
	prev_loss = None
	for a in all_rounds:
		if 'p_value' in a:
			if not prev_loss:
				prev_loss = a['loss']
				continue
			diff = total_p * (prev_loss - a['loss']) / prev_loss
			if diff < threshold:
				return a
	return all_rounds[-1]

def choose_round_diff_read(all_rounds, total_p, read_count, threshold= 1000):
	# max_k = int(total_p * read_count / 1000)
	min_p = threshold/(total_p * read_count)
	# print(min_p)
	# print(max_k)
	rounds = {}
	max_k = 1
	for a in all_rounds:
		max_k = max(max_k, a['k'])
		rounds[a['k']] = a

	new_rounds = []
	prev = None
	for k in range(1,max_k+1):
		a = rounds[k]
		new_rounds.append(a)
		# print("minp: ", min(a['p']))
		if min(a['p']) < min_p:
			return prev
			return new_rounds[:-1]
		# if not prev:
		# 	prev = a
		# 	continue
		# diff = (prev['loss'] - a['loss']) / prev['loss']
		# # print("diff: ", diff)
		# if diff < threshold:
		# 	# return a
		# 	return new_rounds
		prev = a

	return all_rounds[-1]
	# print(max_k)
	return None

def choose_round_loss(all_rounds, total_p, threshold):
	threshold = 1e-5
	for a in all_rounds:
		if 'p_value' in a:
			if a['loss'] * total_p < threshold:
				return a
	return all_rounds[-1]

def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', required=True, help="Input Read File")
	# parser.add_argument('-', '--input', required=True, help="Input Read File")
	parser.add_argument('-f', '--files', nargs='+', required=True, help="file dirs")
	parser.add_argument('-t', '--threshold', required=False, default=1000, help="Threshold")
	parser.add_argument('-o', '--output', required=False, help="Output file")

	args = parser.parse_args()

	placements = {}

	threshold= float(args.threshold)

	with open(args.input, "r") as f:
			lines = f.readlines()
			if len(lines) == 1:
				read_count = int(float(lines[0].split()[0]))
			else:
				# print([float(l.split()[1]) for l in lines])
				read_count = np.sum([float(l.split()[1]) for l in lines])
			# print(read_count)

	# for t in range(16):
	for indir in args.files:
		# print("tree: ", t)
		# indir=os.path.join(args.input, "tree_" + str(t) + "_q1_r2")

		with open(os.path.join(indir, "total_p.txt"), "r") as f:
			total_p = float(f.readlines()[0].split()[0])
		
		if not os.path.exists(os.path.join(indir, "all_rounds.json")):
			continue

		with open(os.path.join(indir, "all_rounds.json")) as f:
			all_rounds = json.load(f)
		fround = None

		if len(all_rounds) > 0:
			fround = choose_round_diff_read(all_rounds, total_p, read_count, threshold = threshold)
		# return
		# fround = None
		# if len(all_rounds) > 0:
		# 	fround = all_rounds[-1]
		# if fround is None:
		# 	print(args.input, t)
		# else:
		# 	# print(len(fround))
		# 	with open(os.path.join(indir, args.output), 'w') as f:
		# 		json.dump(fround, f)


		# print(fround['k'])
		if fround:
			for i in range(len(fround['anchors'])):
				r = fround['p'][i] * total_p * read_count
				if fround['p'][i] < 0:
					r = 0
				if fround['anchors'][i] in placements:
					placements[fround['anchors'][i]] += r
				else:
					placements[fround['anchors'][i]] = r

	# return
	for p in placements:
		print(p, "\t", placements[p])

if __name__ == "__main__":
	main()  