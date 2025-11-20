import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import numpy as np
import seaborn as sns
from upsetplot import UpSet, from_memberships

def compare_runtime(data, no_bowtie2=False, log=True, ymax=None):
    """Compares the runtimes of the benchmark from overview.tsv data

    Args:
        data (pd.DataFrame): overview.tsv loaded as dataframe
    """
    if no_bowtie2:
        df = data[['Bronko_Time', 'BT2_LoFreq_Time', 'BT2_iVar_Time']]
        labels = ['bronko', 'Bowtie2+LoFreq', 'Bowtie2+iVar']
    else:
        df = data[['Bronko_Time','BT2_Time', 'BT2_LoFreq_Time', 'BT2_iVar_Time']]
        labels = ['bronko','Bowtie2 Alone', 'Bowtie2+LoFreq', 'Bowtie2+iVar']
    df_long = df.reset_index().melt(id_vars="index", var_name="Tool", value_name="Runtime")
    plt.figure(figsize=(max(4, len(df.columns)), 5))

    # Boxplot for each column
    # tool_order = ['bronko', 'bowtie2', 'lofreq', 'ivar']
    tool_order = ['Bronko_Time', 'BT2_Time', 'BT2_LoFreq_Time', 'BT2_iVar_Time']
    base_palette = sns.color_palette("colorblind", n_colors=4)
    tmp_palette = [base_palette[1], "gray", base_palette[2], base_palette[0]]
    # tmp_palette = ["seagreen", "cornflowerblue", "moccasin"]
    palette = dict(zip(tool_order, tmp_palette))
    
    ax = sns.barplot(data=df_long, x="Tool", y="Runtime", palette=palette, hue_order=tool_order)
    for container in ax.containers:
        ax.bar_label(container, label_type='edge', fmt="%.2f", fontweight='bold', padding=5)

    # Scatter (strip) plot on top
    sns.stripplot(data=df_long, x="Tool", y="Runtime", color="black", size=4, jitter=True, alpha=0.6)

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.grid(True, linestyle="--", alpha=0.3)
    plt.title("Runtime Comparison")
    plt.ylabel("Runtime (s)")
    ax.set_xlabel("")
    if log:
        plt.yscale("log")
    plt.ylim(bottom=0)
    if ymax:
        plt.ylim(1, ymax)
    ax.set_xticks(ticks=range(len(df.columns)),labels=labels)
    plt.show()
    
def upset_plot_variants(data):
    """Plots upset plots given a major_variants or minor_variants df generated from bench

    Args:
        data (pd.Dataframe): major or minor variants
    """
    df = data
    tool_data = df['tools'].str.split(',')
    tool_data = from_memberships(tool_data)
    fig = plt.figure(figsize=(12, 8))
    upset = UpSet(tool_data, subset_size='count', facecolor="darkblue", shading_color="lightgray", show_counts="%d", sort_by="cardinality").plot(fig=fig)
    plt.suptitle("Variants shared between each method", fontsize=16)
    plt.ylabel("Number of Variants", fontsize=14)
      
    
def upset_plot_stacked(data):
    """Plots UpSet plots of variants called by tools, stacked by TP/FP classification.

    Args:
        data (pd.DataFrame): DataFrame with 'tools' (comma-separated list of tools)
                             and 'truth' (1 = TP, 0 = FP)
    """
    # Split tool memberships
    data = data[['tools', 'truth']]
    memberships = data['tools'].str.split(',')
    
    # Add the stacked variable (must be categorical or string for coloring)
    # Optional: map 1/0 to TP/FP for clearer legend
    data['truth_label'] = data['truth'].map({1: 'True Positive', 0: 'False Positive'})
    
    # Create the multi-index required by UpSet (with truth info for coloring)
    upset_data = from_memberships(memberships, data=data)

    # Initialize the UpSet plot without the default bar chart
    upset = UpSet(
        upset_data,
        subset_size='count',
        facecolor="darkblue", 
        shading_color="lightgray",
        show_counts='%.0f',
        sort_by='cardinality',
        intersection_plot_elements=0
    )

    # Add stacked bars grouped by TP/FP
    upset.add_stacked_bars(
        by="truth_label",                   # column to split bars
        colors={'True Positive': 'Green', 'False Positive': 'Red'},                     # use a colormap (or pass dict)
        title="FP/TP Variant Calls",             # y-axis label
    )

    # Plot
    fig = plt.figure(figsize=(12, 8))
    upset.plot(fig=fig)
    # plt.suptitle("Variants Shared Between Each Method (TP vs FP)")
    plt.show()
    
def compute_stats_with_counts(row, suffix):
    """Adds stats to a row of counts
    """
    bases = ['A', 'C', 'G', 'T']
    counts = [row[f"{b.upper()}_{suffix}"] + row[f"{b.lower()}_{suffix}"] for b in bases]
    depth = sum(counts)
    sorted_counts = sorted(zip(bases, counts), key=lambda x: x[1], reverse=True)
    major, major_count = sorted_counts[0]
    minor, minor_count = sorted_counts[1]
    maf = minor_count / depth if depth > 0 else 0
    return pd.Series([
        depth,
        major,
        major_count,
        minor,
        minor_count,
        maf
    ], index=[
        f'{suffix}_depth',
        f'{suffix}_major',
        f'{suffix}_major_count',
        f'{suffix}_minor',
        f'{suffix}_minor_count',
        f'{suffix}_maf'
    ])

def add_truth_info(minor_counts, truth_loc):
    minor_counts['truth'] = minor_counts['index'].isin(truth_loc).astype(int)
    return minor_counts

def plot_sim_truth(minor_counts_truth, total_TP=1000):
    df = minor_counts_truth[['tools', 'truth']]
    df_exploded = df.assign(tool=df['tools'].str.split(',')).explode('tool')

    counts = df_exploded.groupby(['truth', 'tool']).size().reset_index(name='count')
    counts['truth'] = counts['truth'].map({1: 'True Positives', 0: 'False Positives'})
    
    tool_order = ['bronko', 'lofreq', 'ivar']
    base_palette = sns.color_palette("colorblind", n_colors=4)
    tmp_palette = [base_palette[1], base_palette[2], base_palette[0]]
    # tmp_palette = ["seagreen", "cornflowerblue", "moccasin"]
    palette = dict(zip(tool_order, tmp_palette))
    counts['tool'] = pd.Categorical(counts['tool'], categories=tool_order, ordered=True)
    
    precision_df = counts.pivot(index='tool', columns='truth', values='count').fillna(0)
    precision_df['Precision'] = precision_df['True Positives'] / (
        precision_df['True Positives'] + precision_df['False Positives']
    )   
    precision_df['Recall'] = precision_df['True Positives'] / total_TP
    precision_df['F1'] = 2*((precision_df['Precision']*precision_df['Recall']) / (precision_df['Precision']+precision_df['Recall']))
    precision_df = precision_df.reset_index()
    
    
    # Plot with truth on x-axis, tools as hue (TP and FP groups visually separated)
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 4, figsize=(12, 4), gridspec_kw={'width_ratios': [2, 1, 1, 1]})

    # TP/FP count barplot
    sns.barplot(data=counts, x='truth', y='count', hue='tool', dodge=True, ax=axes[0], palette=palette, hue_order=["bronko", "lofreq", "ivar"])
    axes[0].set_title('FP and TP on simulated dataset', fontsize=16)
    axes[0].set_xlabel('')
    axes[0].set_ylabel('Count', fontsize=14)
    axes[0].legend(title='Tool', fontsize=12)
    for container in axes[0].containers:
        axes[0].bar_label(container, label_type='edge', fontweight='bold')


    # Precision barplot
    sns.barplot(data=precision_df, x='tool', y='Precision', hue='tool', ax=axes[1], order=tool_order, palette=palette, hue_order=["bronko", "lofreq", "ivar"])
    axes[1].set_title('Precision by tool', fontsize=16)
    axes[1].set_ylim(0, 1.05)
    axes[1].set_ylabel('Precision', fontsize=14)
    axes[1].set_xlabel('')
    for container in axes[1].containers:
        axes[1].bar_label(container, label_type='edge', fmt="%.2f", fontweight='bold')
        
    sns.barplot(data=precision_df, x='tool', y='Recall', hue='tool', ax=axes[2], order=tool_order, palette=palette, hue_order=["bronko", "lofreq", "ivar"])
    axes[2].set_title('Recall by tool', fontsize=16)
    axes[2].set_ylim(0, 1.05)
    axes[2].set_ylabel('Recall', fontsize=14)
    axes[2].set_xlabel('')
    for container in axes[2].containers:
        axes[2].bar_label(container, label_type='edge', fmt="%.2f", fontweight='bold')
        
    sns.barplot(data=precision_df, x='tool', y='F1', hue='tool', ax=axes[3], order=tool_order, palette=palette, hue_order=["bronko", "lofreq", "ivar"])
    axes[3].set_title('F1 by tool', fontsize=16)
    axes[3].set_ylim(0, 1.05)
    axes[3].set_ylabel('F1', fontsize=14)
    axes[3].set_xlabel('')
    for container in axes[3].containers:
        axes[3].bar_label(container, label_type='edge', fmt="%.2f", fontweight='bold')

    plt.tight_layout()
    plt.show()
    
def plot_tp_fp_vs_maf(minor_counts_truth):
    maf_thresholds = [0.001, 0.005, 0.01]
    tool_order = ['bronko', 'lofreq', 'ivar']
    base_palette = sns.color_palette("colorblind", n_colors=4)
    tmp_palette = [base_palette[1], base_palette[2], base_palette[0]]
    palette = dict(zip(tool_order, tmp_palette))

    results = []

    # Compute TP/FP counts at each MAF threshold
    for maf in maf_thresholds:
        df = minor_counts_truth[
            (minor_counts_truth['bt2_maf'] > maf) | (minor_counts_truth['bronko_maf'] > maf)
        ][['tools', 'truth']]

        if df.empty:
            continue

        df_exploded = df.assign(tool=df['tools'].str.split(',')).explode('tool')
        counts = df_exploded.groupby(['truth', 'tool']).size().reset_index(name='count')
        counts['truth'] = counts['truth'].map({1: 'TP', 0: 'FP'})
        counts['MAF'] = maf
        results.append(counts)

    all_counts = pd.concat(results)
    all_counts.loc[all_counts['truth'] == "FP", 'count'] *= -1
    print(all_counts)

    # --- Plot ---
    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(6, 4))

    sns.barplot(
        data=all_counts[all_counts['truth'] == 'TP'],
        x='MAF',
        y='count',
        hue='tool',
        palette=palette,
        hue_order=tool_order,
        ax=ax,
        dodge=True
    )
    sns.barplot(
        data=all_counts[all_counts['truth'] == 'FP'],
        x='MAF',
        y='count',
        hue='tool',
        palette=palette,
        hue_order=tool_order,
        ax=ax,
        dodge=True
    )

    # Styling
    ax.axhline(0, color='black', linewidth=1.2)
    ax.set_title('True and False Positives Across MAF Thresholds', fontsize=12, fontweight='bold')
    ax.set_ylabel('Count (TP above / FP below)', fontsize=10, fontweight='bold')
    ax.set_xlabel('Minor Allele Frequency Threshold', fontsize=10, fontweight='bold')
    ax.set_xticklabels([f"{maf}%" for maf in maf_thresholds])
    ax.get_legend().remove()
    ax.tick_params(axis='x', which='major', pad=0, labelsize=10)
    ax.tick_params(axis='y', which='major', pad=0, labelsize=8)

    # Label bars
    for container in ax.containers:
        ax.bar_label(container, fmt="%d", label_type='edge', fontweight='bold', fontsize=8)

    plt.tight_layout()
    plt.show()

def plot_recall_vs_maf(minor_counts_truth):
    maf_thresholds = [0.001, 0.005, 0.01]
    tool_order = ['bronko', 'lofreq', 'ivar']
    base_palette = sns.color_palette("colorblind", n_colors=4)
    tmp_palette = [base_palette[1], base_palette[2], base_palette[0]]
    palette = dict(zip(tool_order, tmp_palette))

    recall_results = []

    # Compute precision at each MAF threshold
    for maf in maf_thresholds:
        df = minor_counts_truth[
            (minor_counts_truth['bt2_maf'] > maf) | (minor_counts_truth['bronko_maf'] > maf) | (minor_counts_truth['truth_maf'] > maf)
        ][['tools', 'truth', 'truth_maf']]
        

        if df.empty:
            continue
        
        total_truths = len(df[~df['truth_maf'].isna()])
        df_exploded = df.assign(tool=df['tools'].str.split(',')).explode('tool')

        # Count true positives per tool
        counts = (
            df_exploded[df_exploded['truth'] == 1]
            .groupby('tool')
            .size()
            .reindex(tool_order, fill_value=0)
            .reset_index(name='TP')
        )

        counts['Recall'] = counts['TP'] / total_truths
        counts['MAF'] = maf
        recall_results.append(counts[['tool', 'Recall', 'MAF']])

    recall_all = pd.concat(recall_results)

    # --- Plot recall grouped barplot ---
    sns.set_theme(style="white")
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.barplot(
        data=recall_all,
        x='MAF',
        y='Recall',
        hue='tool',
        palette=palette,
        hue_order=tool_order,
        ax=ax
    )

    ax.set_title('iSNV Recall', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 1.1)
    ax.set_ylabel('Recall', fontsize=10, fontweight='bold')
    ax.set_xlabel('Minor Allele Frequency Threshold', fontsize=10, fontweight='bold')
    ax.get_legend().remove()
    ax.set_xticklabels([f"{maf*100}%" for maf in maf_thresholds])
    ax.tick_params(axis='x', labelsize=10, pad=0)
    ax.tick_params(axis='y', labelsize=8, pad=0)

    for container in ax.containers:
        ax.bar_label(container, label_type='edge', fmt="%.2f", fontweight='bold', fontsize=8)

    plt.tight_layout()
    plt.show()
    return df

def plot_precision_recall_f1_vs_maf(minor_counts_truth):
    maf_thresholds = [0.001, 0.005, 0.01, 0.015]
    tool_order = ['bronko', 'lofreq', 'ivar']
    base_palette = sns.color_palette("colorblind", n_colors=4)
    tmp_palette = [base_palette[1], base_palette[2], base_palette[0]]
    palette = dict(zip(tool_order, tmp_palette))

    precision_results = []
    recall_results = []

    # --- Compute both metrics across thresholds ---
    for maf in maf_thresholds:
        # -------------------- RECALL --------------------
        # df_recall = minor_counts_truth[
        #     (minor_counts_truth['bt2_maf'] > maf)
        #     | (minor_counts_truth['bronko_maf'] > maf)
        #     | (minor_counts_truth['truth_maf'] > maf)
        # ][['tools', 'truth', 'truth_maf']]
    
        df_recall = minor_counts_truth[
            (minor_counts_truth['truth_maf'] > maf)
        ][['tools', 'truth', 'truth_maf']]

        if df_recall.empty:
            continue

        total_truths = len(df_recall[~df_recall['truth_maf'].isna()])
        df_exploded = df_recall.assign(tool=df_recall['tools'].str.split(',')).explode('tool')

        counts = (
            df_exploded[df_exploded['truth'] == 1]
            .groupby('tool')
            .size()
            .reindex(tool_order, fill_value=0)
            .reset_index(name='TP')
        )

        counts['Recall'] = counts['TP'] / total_truths
        counts['MAF'] = maf
        recall_results.append(counts[['tool', 'Recall', 'MAF']])

        # -------------------- PRECISION --------------------
        df_prec = minor_counts_truth[
            (minor_counts_truth['bt2_maf'] > maf)
            | (minor_counts_truth['bronko_maf'] > maf)
        ][['tools', 'truth']]

        if df_prec.empty:
            continue

        df_exploded = df_prec.assign(tool=df_prec['tools'].str.split(',')).explode('tool')
        counts = df_exploded.groupby(['truth', 'tool']).size().reset_index(name='count')
        counts['truth'] = counts['truth'].map({1: 'True Positives', 0: 'False Positives'})
        counts['tool'] = pd.Categorical(counts['tool'], categories=tool_order, ordered=True)

        precision_df = counts.pivot(index='tool', columns='truth', values='count').fillna(0)
        for col in ['True Positives', 'False Positives']:
            if col not in precision_df.columns:
                precision_df[col] = 0

        precision_df['Precision'] = precision_df.apply(
            lambda x: x['True Positives'] / (x['True Positives'] + x['False Positives'])
            if (x['True Positives'] + x['False Positives']) > 0 else 0,
            axis=1
        )
        precision_df = precision_df.reset_index()
        precision_df['MAF'] = maf
        precision_results.append(precision_df[['tool', 'Precision', 'MAF']])

    precision_all = pd.concat(precision_results)
    recall_all = pd.concat(recall_results)

    # --- Merge to compute F1 ---
    merged = pd.merge(precision_all, recall_all, on=['tool', 'MAF'], how='inner')
    merged['F1'] = 2 * (merged['Precision'] * merged['Recall']) / (
        merged['Precision'] + merged['Recall']
    )

    # --- Plot Precision, Recall, F1 side by side ---
    sns.set_theme(style="white")
    fig, axes = plt.subplots(1, 3, figsize=(12, 3))

    metrics = ['Precision', 'Recall', 'F1']
    titles = ['iSNV Precision', 'iSNV Recall', 'iSNV F1 Score']

    for ax, metric, title in zip(axes, metrics, titles):
        sns.barplot(
            data=merged,
            x='MAF', y=metric, hue='tool',
            palette=palette, hue_order=tool_order, ax=ax
        )
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_ylim(0, 1.1)
        ax.set_ylabel(metric, fontsize=10, fontweight='bold')
        ax.set_xlabel('Minor Allele Frequency Threshold', fontsize=10, fontweight='bold')
        ax.set_xticklabels([f"{maf*100}%" for maf in maf_thresholds])
        ax.tick_params(axis='x', labelsize=10, pad=0)
        ax.tick_params(axis='y', labelsize=8, pad=0)
        ax.get_legend().remove()

        for container in ax.containers:
            ax.bar_label(container, label_type='edge', fmt="%.2f", fontweight='bold', fontsize=8)

    plt.tight_layout()
    plt.show()

    return merged

    
def plot_precision_vs_maf(minor_counts_truth):
    maf_thresholds = [0.001, 0.005, 0.01]
    tool_order = ['bronko', 'lofreq', 'ivar']
    base_palette = sns.color_palette("colorblind", n_colors=4)
    tmp_palette = [base_palette[1], base_palette[2], base_palette[0]]
    palette = dict(zip(tool_order, tmp_palette))

    precision_results = []

    # Compute precision at each MAF threshold
    for maf in maf_thresholds:
        df = minor_counts_truth[
            (minor_counts_truth['bt2_maf'] > maf) | (minor_counts_truth['bronko_maf'] > maf)
        ][['tools', 'truth']]

        if df.empty:
            continue

        df_exploded = df.assign(tool=df['tools'].str.split(',')).explode('tool')
        counts = df_exploded.groupby(['truth', 'tool']).size().reset_index(name='count')
        counts['truth'] = counts['truth'].map({1: 'True Positives', 0: 'False Positives'})
        counts['tool'] = pd.Categorical(counts['tool'], categories=tool_order, ordered=True)

        precision_df = counts.pivot(index='tool', columns='truth', values='count').fillna(0)
        for col in ['True Positives', 'False Positives']:
            if col not in precision_df.columns:
                precision_df[col] = 0

        # Safe precision calculation
        precision_df['Precision'] = precision_df.apply(
            lambda x: x['True Positives'] / (x['True Positives'] + x['False Positives'])
            if (x['True Positives'] + x['False Positives']) > 0 else 0,
            axis=1
        )
        precision_df = precision_df.reset_index()
        precision_df['MAF'] = maf

        precision_results.append(precision_df[['tool', 'Precision', 'MAF']])

    # Combine across thresholds
    precision_all = pd.concat(precision_results)

    # --- Plot combined precision barplot ---
    # --- Precision grouped barplot (styled like your current one) ---
    sns.set_theme(style="white")
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.barplot(
        data=precision_all,
        x='MAF',
        y='Precision',
        hue='tool',
        palette=palette,
        hue_order=tool_order,
        ax=ax
    )

    ax.set_title('iSNV Precision', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 1.1)
    ax.set_ylabel('Precision', fontsize=10, fontweight='bold')
    ax.set_xlabel('Minor Allele Frequency Threshold', fontsize=10, fontweight='bold')
    ax.get_legend().remove()
    ax.set_xticklabels([f"{maf*100}%" for maf in maf_thresholds])
    ax.tick_params(axis='x', which='major', pad=0, labelsize=10)
    ax.tick_params(axis='y', which='major', pad=0, labelsize=8)

    for container in ax.containers:
        ax.bar_label(container, label_type='edge', fmt="%.2f", fontweight='bold', fontsize=8)

    plt.tight_layout()
    plt.show()
    
def plot_tps_vs_maf(minor_counts_truth):
    maf_thresholds = [0.001, 0.005, 0.01]
    tool_order = ['bronko', 'lofreq', 'ivar']
    base_palette = sns.color_palette("colorblind", n_colors=4)
    tmp_palette = [base_palette[1], base_palette[2], base_palette[0]]
    palette = dict(zip(tool_order, tmp_palette))

    results = []

    # Compute TP/FP counts at each MAF threshold
    for maf in maf_thresholds:
        df = minor_counts_truth[
            (minor_counts_truth['bt2_maf'] > maf) | (minor_counts_truth['bronko_maf'] > maf)
        ][['tools', 'truth']]

        if df.empty:
            continue

        df_exploded = df.assign(tool=df['tools'].str.split(',')).explode('tool')
        counts = df_exploded.groupby(['truth', 'tool']).size().reset_index(name='count')
        counts['truth'] = counts['truth'].map({1: 'TP', 0: 'FP'})
        counts['MAF'] = maf
        results.append(counts)

    all_counts = pd.concat(results)
    # all_counts.loc[all_counts['truth'] == "FP", 'count'] *= -1

    # --- Plot combined precision barplot ---
    # --- Precision grouped barplot (styled like your current one) ---
    sns.set_theme(style="white")
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.barplot(
        data=all_counts[all_counts['truth'] == 'TP'],
        x='MAF',
        y='count',
        hue='tool',
        palette=palette,
        hue_order=tool_order,
        ax=ax
    )

    ax.set_title('iSNV True Positives', fontsize=12, fontweight='bold')
    # ax.set_ylim(0, 1.1)
    ax.set_ylabel('Number of TPs', fontsize=10, fontweight='bold')
    ax.set_xlabel('Minor Allele Frequency Threshold', fontsize=10, fontweight='bold')
    ax.get_legend().remove()
    ax.set_xticklabels([f"{maf*100}%" for maf in maf_thresholds])
    ax.tick_params(axis='x', which='major', pad=0, labelsize=10)
    ax.tick_params(axis='y', which='major', pad=0, labelsize=8)
    ax.set_ylim(0, all_counts[all_counts['truth'] == 'TP']['count'].max()*1.1)

    for container in ax.containers:
        ax.bar_label(container, label_type='edge', fmt="%i", fontweight='bold', fontsize=8)

    plt.tight_layout()
    plt.show()
    
def plot_fps_vs_maf(minor_counts_truth):

    maf_thresholds = [0.001, 0.005, 0.01]
    tool_order = ['bronko', 'lofreq', 'ivar']
    base_palette = sns.color_palette("colorblind", n_colors=4)
    tmp_palette = [base_palette[1], base_palette[2], base_palette[0]]
    palette = dict(zip(tool_order, tmp_palette))

    results = []

    for maf in maf_thresholds:
        df = minor_counts_truth[
            (minor_counts_truth['bt2_maf'] > maf) | (minor_counts_truth['bronko_maf'] > maf)
        ][['tools', 'truth']]

        if df.empty:
            continue

        df_exploded = df.assign(tool=df['tools'].str.split(',')).explode('tool')
        counts = df_exploded.groupby(['truth', 'tool']).size().reset_index(name='count')
        counts['truth'] = counts['truth'].map({1: 'TP', 0: 'FP'})
        counts['MAF'] = maf
        results.append(counts)

    all_counts = pd.concat(results, ignore_index=True)

    # --- ensure all tool–MAF combos exist, even with 0 count ---
    full_index = pd.MultiIndex.from_product(
        [maf_thresholds, tool_order, ['FP', 'TP']], names=['MAF', 'tool', 'truth']
    )
    all_counts = (
        all_counts
        .set_index(['MAF', 'tool', 'truth'])
        .reindex(full_index, fill_value=0)
        .reset_index()
    )

    # --- Plot ---
    sns.set_theme(style="white")
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.barplot(
        data=all_counts[all_counts['truth'] == 'FP'],
        x='MAF',
        y='count',
        hue='tool',
        palette=palette,
        hue_order=tool_order,
        ax=ax
    )

    ax.set_title('iSNV False Positives', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of FPs', fontsize=10, fontweight='bold')
    ax.set_xlabel('Minor Allele Frequency Threshold', fontsize=10, fontweight='bold')
    ax.get_legend().remove()
    ax.set_xticklabels([f"{maf*100}%" for maf in maf_thresholds])
    ax.tick_params(axis='x', which='major', pad=0, labelsize=10)
    ax.tick_params(axis='y', which='major', pad=0, labelsize=8)
    ax.set_ylim(0, all_counts[all_counts['truth'] == 'FP']['count'].max() * 1.1)

    for container in ax.containers:
        ax.bar_label(container, label_type='edge', fmt="%i", fontweight='bold', fontsize=8)

    plt.tight_layout()
    plt.show()

# Compute IV stats
def get_columns(variants):
    """Appends a bunch of useful columns to the major_variants.tsv or minor_variants.tsv dataframe
        including the major/minor alleles, their frequencies, etc.
    """
    variants = variants.join(variants.apply(compute_stats_with_counts, axis=1, suffix='bronko'))
    variants = variants.join(variants.apply(compute_stats_with_counts, axis=1, suffix='bt2'))

    # Compute differences
    variants['Depth Difference'] = variants['bt2_depth'] - variants['bronko_depth']
    variants['maf difference'] = variants['bt2_maf'] - variants['bronko_maf']

    # Reorder columns (optional)
    desired_order = [
        'Sample', 'reference', 'index', 'ref', 'alt', 'tools',
        'Depth Difference', 'maf difference',
        'bronko_depth', 'bronko_major', 'bronko_major_count','bronko_minor', 'bronko_minor_count','bronko_maf', 'A_bronko', 'C_bronko', 'G_bronko', 'T_bronko', 'a_bronko', 'c_bronko', 'g_bronko', 't_bronko',
        'bt2_depth', 'bt2_major', 'bt2_major_count','bt2_minor', 'bt2_minor_count', 'bt2_maf', 'A_bt2', 'C_bt2', 'G_bt2', 'T_bt2', 'a_bt2', 'c_bt2', 'g_bt2', 't_bt2'
    ]

    variants = variants[desired_order]
    return variants


def plot_maf_comparison(minor_all_stats, no_truth=False):
    """
    Plot two side-by-side MAF comparison scatterplots:
    - Left: full range (0–0.5)
    - Right: zoomed in (0–0.1)

    Parameters:
        minor_all_stats (pd.DataFrame): DataFrame with columns 'bt2_maf', 'bronko_maf', and 'tools'
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    minor_all_stats = minor_all_stats.sort_values(by='tools')
    if not no_truth:
        minor_all_stats['truth'] = minor_all_stats['truth'].map({1: 'TP', 0: 'FP'})


    # Full-range plot
    if no_truth:
        sns.scatterplot(
            data=minor_all_stats,
            x='bt2_maf', y='bronko_maf',
            hue='tools', ax=axes[0],
            palette="colorblind", s=50, edgecolor="black", linewidth=0.5
        )
    else:
        sns.scatterplot(
            data=minor_all_stats,
            x='bt2_maf', y='bronko_maf',
            hue='tools', style='truth', style_order=['TP', 'FP'],
            ax=axes[0], palette="colorblind", s=50, edgecolor="black", linewidth=0.5
        )

    # Reference lines
    axes[0].plot([0, 1], [0, 1], 'k--', alpha=0.6, zorder=0)  # diagonal
    axes[0].axhline(0, color='gray', lw=0.8, alpha=0.7, zorder=0)
    axes[0].axvline(0, color='gray', lw=0.8, alpha=0.7, zorder=0)

    # Limits
    axes[0].set_xlim([0, 0.5])
    axes[0].set_ylim([0, 0.5])

    # Labels + title
    axes[0].set_xlabel("Bowtie2 Pileup Minor Allele Frequency", fontsize=14, fontweight="bold")
    axes[0].set_ylabel("bronko Minor Allele Frequency", fontsize=14, fontweight="bold")
    axes[0].set_title("MAF Comparison (0–0.5)", fontsize=16, fontweight="bold")

    # Grid styling
    axes[0].grid(True, linestyle="--", alpha=0.5)

    # Legend (move outside, only if you want to keep it)
    axes[0].get_legend().remove()


     # Zoomed-in plot
    if no_truth:
        sns.scatterplot(
            data=minor_all_stats,
            x='bt2_maf', y='bronko_maf',
            hue='tools', ax=axes[1],
            palette="colorblind", s=50, edgecolor="black", linewidth=0.5
        )
    else:
        sns.scatterplot(
            data=minor_all_stats,
            x='bt2_maf', y='bronko_maf',
            hue='tools', style='truth', style_order=['TP', 'FP'],
            ax=axes[1], palette="colorblind", s=50, edgecolor="black", linewidth=0.5
        )

    # Get zoomed-in axis range
    max_count = max(
        minor_all_stats['bt2_maf'].quantile(0.95),
        minor_all_stats['bronko_maf'].quantile(0.95), 
        0.03
    )

    # Reference lines
    axes[1].plot([0, 1], [0, 1], 'k--', alpha=0.6, zorder=0)
    axes[1].axhline(0, color='black', linewidth=1)
    axes[1].axvline(0, color='black', linewidth=1)

    # Axis limits
    axes[1].set_xlim([0, max_count])
    axes[1].set_ylim([0, max_count])

    # Labels and title
    axes[1].set_xlabel("Bowtie2 Pileup Minor Allele Frequency", fontsize=12, weight="bold")
    axes[1].set_ylabel("", fontsize=12)
    axes[1].set_title(f"MAF Comparison (0–{round(max_count, 2)})", fontsize=14, weight="bold")

    # Ticks
    axes[1].tick_params(axis='both', which='major', labelsize=10)

    # Grid for readability
    axes[1].grid(True, linestyle="--", alpha=0.6)

    # Remove duplicate legend if shared with first subplot
    axes[1].get_legend().remove()
    # Get handles/labels for both legends
    handles_tools, labels_tools = axes[0].get_legend_handles_labels()

    # If you used style='truth', seaborn merges legends into one —
    # so we split them manually by checking unique labels
    unique_labels = list(dict.fromkeys(labels_tools))  # preserve order
    unique_handles = []
    for lbl in unique_labels:
        idx = labels_tools.index(lbl)
        unique_handles.append(handles_tools[idx])

    # Separate "Tools" and "Truth"
    tools_labels = [lbl for lbl in unique_labels if lbl not in ["TP", "FP"]]
    tools_handles = [unique_handles[unique_labels.index(lbl)] for lbl in tools_labels]
    if tools_labels[0] == 'tools':
        tools_labels = tools_labels[1:]
        tools_handles = tools_handles[1:]
        
    if tools_labels[-1] == 'truth':
        tools_labels = tools_labels[:-1]
        tools_handles = tools_handles[:-1]

    truth_labels = [lbl for lbl in unique_labels if lbl in ["TP", "FP"]]
    truth_handles = [unique_handles[unique_labels.index(lbl)] for lbl in truth_labels]

    # Remove legends from axes
    for ax in axes:
        if ax.get_legend() is not None:
            ax.get_legend().remove()

    # Add figure-level legends
    legend1 = fig.legend(
        tools_handles, tools_labels,
        title="Tools",
        loc="upper left",
        bbox_to_anchor=(1.0, 0.76),
        borderaxespad=0
    )

    if not no_truth:
        legend2 = fig.legend(
            truth_handles, truth_labels,
            title="Truth",
            loc="upper left",
            bbox_to_anchor=(0.99, 0.43),  # stack underneath Tools legend
            borderaxespad=0
        )
        plt.setp(legend2.get_title(), fontweight="bold")


    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.setp(legend1.get_title(), fontweight="bold")

    plt.tight_layout()
    plt.show()
    
    
def plot_count_comparison(minor_all_stats, adjust=False, k=None, b=None, e=0.01, per=0.1, no_truth=False):
    """
    Plot two side-by-side MAF comparison scatterplots:
    - Left: full range (0–0.5)
    - Right: zoomed in (0–0.1)

    Parameters:
        minor_all_stats (pd.DataFrame): DataFrame with columns 'bt2_maf', 'bronko_maf', and 'tools'
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    
    if not no_truth:
        minor_all_stats = minor_all_stats.sort_values(by='tools')
        minor_all_stats['truth'] = minor_all_stats['truth'].map({1: 'TP', 0: 'FP'})


    if adjust and k and b:
        minor_all_stats['bronko_minor_count'] = (minor_all_stats['bronko_minor_count'] / (100 - k + 2*b + per) * 100).astype(int)
        minor_all_stats['bronko_major_count'] = (minor_all_stats['bronko_major_count'] / (100 - k + 2*b + e) * 100).astype(int)

    # Full-range plot
    if no_truth:
        sns.scatterplot(
            data=minor_all_stats,
            x='bt2_minor_count', y='bronko_minor_count',
            hue='tools', ax=axes[0][0]
        )
    else:
        sns.scatterplot(
            data=minor_all_stats,
            x='bt2_minor_count', y='bronko_minor_count',
            hue='tools', ax=axes[0][0], style='truth', style_order=['TP', 'FP']
        )
    max_count = max(minor_all_stats['bt2_minor_count'].max(), minor_all_stats['bronko_minor_count'].max())
    axes[0][0].plot([0, max_count], [0, max_count], 'k-', alpha=0.75, zorder=0)
    axes[0][0].axhline(0, color='black')
    axes[0][0].axvline(0, color='black')
    axes[0][0].set_xlim([0, max_count])
    axes[0][0].set_ylim([0, max_count])
    axes[0][0].set_ylabel("bronko Pileup Minor Allele Count")
    axes[0][0].set_xlabel("Bowtie2 Pileup Minor Allele Count")
    axes[0][0].set_title("Minor Allele Count Comparison")
    
        # Full-range plot
    if no_truth:
        sns.scatterplot(
            data=minor_all_stats,
            x='bt2_minor_count', y='bronko_minor_count',
            hue='tools', ax=axes[0][1]
        )
    else: 
        sns.scatterplot(
            data=minor_all_stats,
            x='bt2_minor_count', y='bronko_minor_count',
            hue='tools', ax=axes[0][1], style='truth', style_order=['TP', 'FP']
        )
    
    max_count = max(minor_all_stats['bt2_minor_count'].quantile(q=0.95).max(), minor_all_stats['bronko_minor_count'].quantile(q=0.95).max())
    axes[0][1].plot([0, max_count], [0, max_count], 'k-', alpha=0.75, zorder=0)
    axes[0][1].axhline(0, color='black')
    axes[0][1].axvline(0, color='black')
    axes[0][1].set_xlim([0, max_count])
    axes[0][1].set_ylim([0, max_count])
    axes[0][1].set_ylabel("bronko Pileup Minor Allele Count")
    axes[0][1].set_xlabel("Bowtie2 Pileup Minor Allele Count")
    axes[0][1].set_title("Minor Allele Count Comparison")
    
    # Zoomed-in plot
    if no_truth:
        sns.scatterplot(
            data=minor_all_stats,
            x='bt2_major_count', y='bronko_major_count',
            hue='tools', ax=axes[1][0], legend=False
        )
    else:
        sns.scatterplot(
            data=minor_all_stats,
            x='bt2_major_count', y='bronko_major_count',
            hue='tools', ax=axes[1][0], legend=False, style='truth', style_order=['TP', 'FP']
        )
    max_count = max(minor_all_stats['bt2_major_count'].max(), minor_all_stats['bronko_major_count'].max())
    axes[1][0].plot([0, max_count], [0, max_count], 'k-', alpha=0.75, zorder=0)
    axes[1][0].axhline(0, color='black')
    axes[1][0].axvline(0, color='black')
    axes[1][0].set_xlim([0, max_count])
    axes[1][0].set_ylim([0, max_count])
    axes[1][0].set_xlabel("Bowtie2 Pileup Major Allele Count")
    axes[1][0].set_ylabel("bronko Pileup Minor Allele Count")
    axes[1][0].set_title("Major Allele Count Comparison")
    
    
            # Full-range plot
    if no_truth:
        sns.scatterplot(
            data=minor_all_stats,
            x='bt2_major_count', y='bronko_major_count',
            hue='tools', ax=axes[1][1]
        )    
    else:
        sns.scatterplot(
            data=minor_all_stats,
            x='bt2_major_count', y='bronko_major_count',
            hue='tools', ax=axes[1][1], style='truth', style_order=['TP', 'FP']
        )
        
    
    max_count = max(minor_all_stats['bt2_major_count'].quantile(q=0.95).max(), minor_all_stats['bronko_major_count'].quantile(q=0.95).max())
    axes[1][1].plot([0, max_count], [0, max_count], 'k-', alpha=0.75, zorder=0)
    axes[1][1].axhline(0, color='black')
    axes[1][1].axvline(0, color='black')
    axes[1][1].set_xlim([0, max_count])
    axes[1][1].set_ylim([0, max_count])
    axes[1][1].set_ylabel("bronko Pileup Major Allele Count")
    axes[1][1].set_xlabel("Bowtie2 Pileup Major Allele Count")
    axes[1][1].set_title("Major Allele Count Comparison")

    # Shared legend above plots
    # handles, labels = axes[0].get_legend_handles_labels()
    # fig.legend(handles, labels, title="Variant Identified", loc='right', bbox_to_anchor=(1.12, 0.5))

    plt.tight_layout()
    plt.show()
    
