#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path


def parse_exp_id(exp_id):
    dataset_id = exp_id.split('__cells-')[0]
    cells = re.search(r'__cells-([^_]+)__', exp_id)
    rep = re.search(r'__rep-(\d+)__', exp_id)
    return dataset_id, (cells.group(1) if cells else 'all'), int(rep.group(1)) if rep else 1


def read_old_seed(checkpoint_root, old_exp):
    metrics = checkpoint_root.parent / f'{old_exp}.metrics.tsv'
    if not metrics.exists():
        return ''
    with metrics.open(newline='') as fh:
        reader = csv.DictReader(fh, delimiter='	')
        for row in reader:
            seed = row.get('seed', '')
            if seed:
                return seed
    return ''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--checkpoint-root', default='/home/dzf/cellchat_acceleration/results/runs/checkpoints')
    ap.add_argument('--out', default='results/sparse_exact_checkpoint_grid.csv')
    ap.add_argument('--algorithm', default='sparse_exact')
    ap.add_argument('--limit', type=int, default=0)
    args = ap.parse_args()

    root = Path(args.checkpoint_root)
    paths = sorted(root.glob('*__engine-both__ablation-full/prepared_cellchat.rds'))
    rows = []
    for i, path in enumerate(paths, start=1):
        old_exp = path.parent.name
        dataset_id, cells, rep = parse_exp_id(old_exp)
        seed = read_old_seed(root, old_exp) or str(100000 + i)
        # Keep the old experiment id visible but avoid overwriting existing dense checkpoints.
        exp_id = f'{old_exp}__algorithm-{args.algorithm}'
        rows.append({
            'experiment_id': exp_id,
            'dataset_id': dataset_id,
            'input_path': '',
            'prepared_cellchat': str(path),
            'n_cells': cells,
            'repeat': rep,
            'seed': seed,
            'engine': 'accelerated',
            'ablation': 'full',
            'accel_algorithm': args.algorithm,
            'label_col': 'auto',
        })
        if args.limit and len(rows) >= args.limit:
            break

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w', newline='') as fh:
        fieldnames = ['experiment_id', 'dataset_id', 'input_path', 'prepared_cellchat',
                      'n_cells', 'repeat', 'seed', 'engine', 'ablation',
                      'accel_algorithm', 'label_col']
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f'checkpoint_rows={len(rows)} out={out}')


if __name__ == '__main__':
    main()
