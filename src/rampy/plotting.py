# -*- coding: utf-8 -*-
#############################################################################
#Copyright (c) 2018-2025 Charles Le Losq
#
# Licence GNU-GPL
#
#
#############################################################################
import plotly.graph_objects as go
import numpy as np

def plot_spectrum(
    x,
    y,
    fig=None,
    baselines=None,
    baseline_labels=None,
    smoothed_signals=None,
    smoothed_labels=None,
    label="Spectrum",
    color=None,
    roi=None,
    xaxis_title=None,
    yaxis_title=None,
    title=None
):
    """Plots a spectrum with optional baselines, smoothed signals, and regions of interest (bir) using Plotly.

    This function adds a single spectrum to a Plotly figure. It supports overlaying baselines, smoothed signals, 
    and highlighting regions of interest (bir). If no figure is provided, a new Plotly figure is created.

    Args:
        x (array-like): The x-axis values for the spectrum.
        y (array-like): The y-axis values for the spectrum.
        fig (plotly.graph_objects.Figure, optional): The Plotly figure to which the spectrum will be added. 
            If None, a new figure will be created. Default is None.
        baselines (list of array-like, optional): List of baseline arrays to overlay on the spectrum. Default is None.
        baseline_labels (list of str, optional): Labels for each baseline. Default is None.
        smoothed_signals (list of array-like, optional): List of smoothed signal arrays to overlay on the spectrum. Default is None.
        smoothed_labels (list of str, optional): Labels for each smoothed signal. Default is None.
        label (str, optional): Label for the spectrum. Default is "Spectrum".
        color (str or tuple, optional): Base color for the spectrum. If None, a default color will be chosen. Default is None.
        bir (list or ndarray, optional): Regions of interest to highlight on the plot. Should be a list or array 
            of shape `(n_regions, 2)`, where each row specifies the start and end of a region 
            (e.g., `[[x_start1, x_end1], [x_start2, x_end2]]`). Default is None.
        xaxis_title (str, optional): Title for the x-axis. If None, no specific title is set. Default is None.
        yaxis_title (str, optional): Title for the y-axis. If None, no specific title is set. Default is None.
        title (str, optional): Title for the entire plot. If None, no title is set. Default is None.

    Returns:
        plotly.graph_objects.Figure: The updated Plotly figure object containing the plotted spectrum.

    Examples:
        Plot a simple spectrum:

        >>> import numpy as np
        >>> import plotly.graph_objects as go
        >>> x = np.linspace(0, 10, 500)
        >>> y = np.sin(x) + np.random.normal(0, 0.1, len(x))
        >>> fig = plot_spectrum(
                x=x,
                y=y,
                label="Sample Spectrum",
                xaxis_title="Wavenumber (cm⁻¹)",
                yaxis_title="Intensity (a.u.)",
                title="Sample Spectrum Analysis"
            )
        >>> fig.show()

        Plot a spectrum with baseline and smoothed signal:

        >>> baseline = 0.5 * np.ones_like(x)
        >>> smoothed = np.convolve(y - baseline, np.ones(10)/10, mode='same')
        >>> fig = plot_spectrum(
        >>>        x=x,
        >>>        y=y,
        >>>        baselines=[baseline],
        >>>        baseline_labels=["Constant Baseline"],
        >>>        smoothed_signals=[smoothed],
        >>>        smoothed_labels=["Moving Average"],
        >>>        xaxis_title="Wavenumber (cm⁻¹)",
        >>>        yaxis_title="Intensity (a.u.)"
        >>>    )
        >>> fig.show()
    
    Notes:
        - All input arrays are automatically converted to vectors using ``np.asarray(...).flatten()`` to ensure compatibility with Plotly.
        - Baselines are displayed as dashed lines with distinct colors.
        - Smoothed signals are displayed as dotted lines with distinct colors.
        - Regions of interest (BIR) are highlighted with semi-transparent yellow rectangles.
        - The function returns the figure object, allowing for further customization if needed.
    """
    
    if fig is None:  # Create the figure if it is not provided
        fig = go.Figure()

    # Ensure x and y are vectors
    x = np.asarray(x).flatten()
    y = np.asarray(y).flatten()

    # Add the original spectrum
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="lines",
            name=f"{label} - Original",
            line=dict(color=color or "blue", width=2)
        )
    )
    
    # Add baselines if provided
    if baselines is not None:
        if baseline_labels is None:
            baseline_labels = [f"{label} - Baseline {i+1}" for i in range(len(baselines))]
        
        for i, baseline in enumerate(baselines):
            # Ensure baseline is a vector
            baseline = np.asarray(baseline).flatten()
            
            # Use distinct colors for each baseline
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=baseline,
                    mode="lines",
                    name=baseline_labels[i],
                    line=dict(dash="dash", color=f"rgba({(i*50)%255}, {(i*100)%255}, {(i*150)%255}, 0.8)", width=1.5)
                )
            )
    
    # Add smoothed signals if provided
    if smoothed_signals is not None:
        if smoothed_labels is None:
            smoothed_labels = [f"{label} - Smoothed {i+1}" for i in range(len(smoothed_signals))]
        
        for i, smoothed_signal in enumerate(smoothed_signals):
            # Ensure smoothed signal is a vector
            smoothed_signal = np.asarray(smoothed_signal).flatten()
            
            # Use distinct colors for each smoother
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=smoothed_signal,
                    mode="lines",
                    name=smoothed_labels[i],
                    line=dict(dash="dot", color=f"rgba({(i*100+50)%255}, {(i*150+50)%255}, {(i*200+50)%255}, 0.8)", width=1.5)
                )
            )
    
    # Highlight regions of interest (roi) if provided
    if roi is not None:
        y_min = min(y)
        y_max = max(y)
        y_range = y_max - y_min
        
        for i, region in enumerate(roi):
            fig.add_shape(
                type="rect",
                x0=region[0],
                x1=region[1],
                y0=y_min - 0.05 * y_range,
                y1=y_max + 0.05 * y_range,
                fillcolor="yellow",
                opacity=0.2,
                layer="below",
                line_width=0
            )
            # Add annotation for each roi region
            fig.add_annotation(
                x=np.mean(region),
                y=y_max + 0.02 * y_range,
                text=f"ROI {i+1}",
                showarrow=False,
                font=dict(size=10),
                align="center"
            )
    
    # Update layout with axis titles and plot title
    layout_updates = {}
    
    if xaxis_title is not None:
        layout_updates["xaxis"] = dict(title=dict(text=xaxis_title))
    
    if yaxis_title is not None:
        layout_updates["yaxis"] = dict(title=dict(text=yaxis_title))
    
    if title is not None:
        layout_updates["title"] = title
    
    if layout_updates:
        fig.update_layout(**layout_updates)
    
    return fig
