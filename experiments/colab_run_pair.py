"""Head-to-head on one machine: KaPoCE and memetic_sa with the same budget."""
import json, os, sys, time, threading
sys.path.insert(0, '/content/cc'); sys.path.insert(0, '/content/cc/experiments')
os.environ['KAPOCE_BIN'] = '/content/kapoce/build/ClusterEditing'
import numpy as np
import datasets as D, baselines as B, ccbench as cc
from ccbench.memetic import memetic_sa, memetic_twin

def main(graphs, T, seed, tag):
    os.makedirs('/content/results', exist_ok=True)
    for name in graphs:
        g = D.load(name)
        res = {'graph': name, 'n': g.n, 'm': g.m, 'T': T, 'seed': seed, 'tag': tag}
        box = {}
        def kp():
            lab, t = B.kapoce(g, T)
            box['kapoce'] = (cc.cost(g, lab), t)
        th = threading.Thread(target=kp); th.start()
        t0 = time.time(); hist = []
        algo = memetic_sa if tag.startswith('e1') else memetic_twin
        lab, c = algo(g, T, rng=seed, history=hist)
        res['ours'] = int(cc.cost(g, lab)); res['ours_time'] = time.time() - t0
        th.join()
        res['kapoce'], res['kapoce_time'] = int(box['kapoce'][0]), box['kapoce'][1]
        res['hist'] = hist[::max(1, len(hist) // 50)]
        with open(f'/content/results/{tag}_{name}_{seed}.json', 'w') as fh:
            json.dump(res, fh)
        print(res['graph'], 'ours', res['ours'], 'kapoce', res['kapoce'], flush=True)

if __name__ == '__main__':
    main(sys.argv[4:], float(sys.argv[1]), int(sys.argv[2]), sys.argv[3])
