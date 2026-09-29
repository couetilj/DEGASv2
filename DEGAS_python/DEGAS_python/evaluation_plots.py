"""Optional Matplotlib evaluation figures; no plotting imports during training."""
from pathlib import Path
import numpy as np
from .validation import METRICS

COLORS = ('#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#000000')


def plot_validation(metrics, selection, *, title='DEGAS validation', output_dir=None):
    """Return AUROC and five-metric figures, optionally save PNG/PDF and CSV.

    Takes the two tables from validation.select_one_se. Study curves are faint,
    means are equally study weighted. The vertical marker is exactly one SE,
    not a confidence interval. Training/tuning results are labelled as such.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter, NullLocator
    if metrics.duplicated(['study', 'size']).any() or selection['size'].duplicated().any():
        raise ValueError('Duplicate study/size or selection rows')
    sizes = sorted(selection['size'])
    if not sizes or min(sizes) <= 0:
        raise ValueError('Need positive feature sizes')
    if any(sorted(p['size']) != sizes for _, p in metrics.groupby('study')):
        raise ValueError('Every study must cover every selection size')
    s = selection.sort_values('size')
    means = metrics.groupby('size')[list(METRICS)].mean().loc[sizes]
    if not np.isfinite(metrics[list(METRICS)].to_numpy()).all():
        raise ValueError('Metrics must be finite')
    best = s.iloc[int(np.argmax(s.mean_AUROC.to_numpy()))]
    cutoff = float(best.mean_AUROC - best.SE)
    if (not np.isfinite(s[['mean_AUROC', 'SE', 'cutoff']].to_numpy()).all()
            or (s.SE < 0).any() or not np.allclose(s.cutoff, cutoff)
            or not np.allclose(means.AUROC, s.mean_AUROC)):
        raise ValueError('Selection mean/SE/cutoff do not match study metrics')
    with plt.rc_context({'font.size': 12, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42}):
        fig = plt.figure(figsize=(10, 7.5), layout='constrained')
        grid = fig.add_gridspec(3, 1, height_ratios=[8, 1.4, .8])
        ax = fig.add_subplot(grid[0])
        legend_axis = fig.add_subplot(grid[1])
        caption_axis = fig.add_subplot(grid[2])
        legend_axis.axis('off')
        caption_axis.axis('off')
        panels, axes = plt.subplots(2, 3, figsize=(15, 9), layout='constrained')
        for i, (study, part) in enumerate(metrics.groupby('study', sort=True)):
            part = part.sort_values('size')
            color = COLORS[i % len(COLORS)]
            ax.plot(part['size'], part.AUROC, 'o-', color=color, alpha=.35, label=str(study))
            for axis, metric in zip(axes.flat, METRICS):
                axis.plot(part['size'], part[metric], 'o-', color=color, alpha=.5, label=str(study))
        ax.plot(s['size'], s.mean_AUROC, 'o-', color='#111111', lw=2.2, label='Equal-study mean')
        ax.scatter([best['size']], [best.mean_AUROC], s=180, facecolors='none',
                   edgecolors='#E69F00', linewidths=2, zorder=5, label=f'Best mean ({int(best["size"])} genes)')
        chosen = s[s.selected]
        ax.scatter(chosen['size'], chosen.mean_AUROC, marker='s', s=75,
                   facecolors='none', edgecolors='#009E73', label='Within 1 SE')
        ax.axhline(.5, color='.6', ls=':', label='Chance AUROC')
        ax.axhline(cutoff, color='.2', ls='--', label=f'One-SE cutoff ({cutoff:.4f})')
        ax.errorbar(best['size'], best.mean_AUROC, yerr=[[best.SE], [0]], color='.15', capsize=5)
        ax.annotate('1 SE', (best['size'], (best.mean_AUROC+cutoff)/2),
                    xytext=(-12 if best['size'] == max(sizes) else 12, -12),
                    ha='right' if best['size'] == max(sizes) else 'left', textcoords='offset points')
        ax.set(title=title, ylabel='Held-out patient AUROC', ylim=(max(0, min(.5, metrics.AUROC.min(), cutoff) - .08), 1.04))
        legend_axis.legend(*ax.get_legend_handles_labels(), loc='center', ncol=3, frameon=False, fontsize=10)
        labels = ('AUROC ↑', 'Average precision ↑', 'Log loss ↓', 'Brier score ↓', 'Balanced accuracy ↑')
        for axis, metric, label in zip(axes.flat, METRICS, labels):
            axis.plot(means.index, means[metric], 'o-', color='#111111', lw=2, label='Equal-study mean')
            axis.set(ylabel=label)
            if metric in ('AUROC', 'balanced_accuracy'):
                axis.axhline(.5, color='.6', ls=':')
        axes.flat[-1].axis('off')
        handles, names = axes.flat[0].get_legend_handles_labels()
        axes.flat[-1].legend(handles, names, loc='center', frameon=False)
        for axis in [ax, *list(axes.flat)[:-1]]:
            axis.set_xscale('log')
            axis.set_xlabel('Feature-set size')
            ticks = sizes
            if len(sizes) > 8:
                # Keep endpoints and thin crowded labels; every candidate is plotted.
                gap = (np.log10(sizes[-1])-np.log10(sizes[0])) / (10 if axis is ax else 6)
                ticks = [sizes[0]]
                for size in sizes[1:-1]:
                    if np.log10(size/ticks[-1]) >= gap:
                        ticks.append(size)
                if len(ticks) > 1 and np.log10(sizes[-1]/ticks[-1]) < gap:
                    ticks.pop()
                ticks.append(sizes[-1])
            axis.set_xticks(ticks)
            axis.xaxis.set_major_formatter(ScalarFormatter())
            axis.xaxis.set_minor_locator(NullLocator())
            axis.tick_params(axis='x', rotation=45 if len(sizes) > 8 else 35,
                             labelsize=9 if len(sizes) > 8 else 12)
        panels.suptitle(title + ' — tuning metrics')
        caption_axis.text(.5, .5, 'Paired patient bootstrap within study × diagnosis; SE is not a 95% CI.\n'
                          'Size selection uses validation data; assess final performance on the untouched holdout.',
                          ha='center', va='center', fontsize=10, transform=caption_axis.transAxes)
    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        for name, figure in [('auroc_selection', fig), ('evaluation_metrics', panels)]:
            for ext in ('png', 'pdf'):
                figure.savefig(out / f'{name}.{ext}', dpi=200, bbox_inches='tight')
        metrics.to_csv(out / 'metrics_by_study.csv', index=False)
        selection.to_csv(out / 'feature_size_selection.csv', index=False)
    return fig, panels
