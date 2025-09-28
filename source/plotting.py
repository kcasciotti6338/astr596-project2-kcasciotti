# plotting.py

'''
Plotting functions

Orbits, ___ vs time, ____ vs N, evol. screenshote

'''

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from constants import RSUN
from matplotlib.animation import FuncAnimation, PillowWriter

def style_plot(ax=None, figsize=(8, 6), dpi=100):
    """
    Create figure with consistent style settings
    
    Parameters
    ----------
    ax : matplotlib axis, optional
        If provided, plot on this axis. If None, create new figure
    figsize : tuple
        Figure size in inches
    dpi : int
        Resolution for display
    
    Returns
    -------
    fig, ax : matplotlib objects
        Figure and single axis for plotting
    
    Examples
    --------
    >>> fig, ax = setup_plot()
    >>> ax.plot(x, y)
    >>> plt.show()
    """
    
    if ax == None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    else:
        fig = ax.figure
        
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title('title', fontsize=16)
    ax.tick_params(which="both", labelsize=8)
    fig.tight_layout()
    ax.grid(True, alpha=0.5, which="both")
    
    return fig, ax       

def plot_orbit(results, labl, titl, masses=None, ax=None, **kwargs):
    """
    Create a properly formatted graph of the orbit
    
    Parameters
    ----------
    results : dataclass
        - ts : times
        - ps : positions
    labl : label for legend
    titl : label for title
    masses : array of masses (optional)
    ax : matplotlib axis, optional
        If provided, plot on this axis. If None, create new figure
    **kwargs : dict
        Additional plotting parameters. Common options:
        - cmap : str, colormap name (default: 'viridis')
        - s : float, marker size (default: 20)
        - alpha : float, transparency (default: 0.7)
        - figsize : tuple, figure dimensions (only used if ax is None)
        
    Returns
    -------
    fig, ax : matplotlib objects
        Figure and axis used for plotting
        
    """
    
    
    N, d = results.ps.shape[1:3]
    
    if ax is None:
        figsize=kwargs.pop('figsize', (8, 8))
        if d == 2:
            fig, ax = plt.subplots(figsize=figsize)
        else: 
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111, projection='3d')
    else:
        fig = ax.get_figure()
    
    cmap = kwargs.get('cmap', 'coolwarm_r')
    alpha = kwargs.get('alpha', 0.7)
    
    if masses is None:
        s = kwargs.get('s', 20)
    else:
        norm = (masses - np.min(masses)) / (np.max(masses) - np.min(masses)) 
        s = 10 + norm*200
    
    base = plt.get_cmap('tab10')
    palette = [base(i % 10) for i in range(N)]
    body_colors = {j: palette[j] for j in range(N)}
    
    for j in range(N):
        col = body_colors[j]
        if d == 2:
            x, y = results.ps[:, j, 0], results.ps[:, j, 1]
            ax.plot(x, y, label=f'{labl} - {j}', **kwargs)
        else:
            x, y, z = results.ps[:, j, 0], results.ps[:, j, 1], results.ps[:, j, 2]
            ax.plot(x, y, z, label=labl, **kwargs)
    
    cols = [body_colors[j] for j in range(N)]    
    if d == 2:    
        ax.scatter(results.ps[-1, :, 0], results.ps[-1, :, 1], **kwargs) 
    else:
        ax.scatter(results.ps[-1, :, 0], results.ps[-1, :, 1], results.ps[-1, :, 2], **kwargs) 
       
    fig, ax = style_plot(ax)
    ax.set_xlabel('X Position (AU)')
    ax.set_ylabel('Y Position (AU)')
    ax.set_title(titl)
    
    return fig, ax 

def plot_snapshot(system, a=100, m=None, titl=None, ax=None, **kwargs):
    """
    Create a properly formatted graph of a snapshot of the evolution
    
    Parameters
    ----------
    system : N-body object
        - position
        - mass
        - dimensions
    a : half-mass radius (AU), optional
    titl : title, optional
    ax : matplotlib axis, optional
        If provided, plot on this axis. If None, create new figure
    **kwargs : dict
        Additional plotting parameters. Common options:
        - cmap : str, colormap name (default: 'viridis')
        - s : float, marker size (default: 20)
        - alpha : float, transparency (default: 0.7)
        - figsize : tuple, figure dimensions (only used if ax is None)
        
    Returns
    -------
    fig, ax : matplotlib objects
        Figure and axis used for plotting
        
    """
    
    if ax is None:
        figsize = kwargs.pop('figsize', (8, 6))
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
    
    cmap = kwargs.get('cmap', 'coolwarm_r')
    s = kwargs.get('s', 20)
    alpha = kwargs.get('alpha', 0.7)
    
    if m is None: 
        m = np.asarray(system.ms, dtype=float).ravel()
    diff = np.ptp(m)
    if diff > 0:
        norm = (m - m.min()) / diff
    else:
        norm = np.zeros_like(m)           

    sM = 2 + norm * 100                   
    sizes = kwargs.pop('s', sM)
    
    Rmax = 2*a
    
    if system.dimensions == 2:
        x, y = system.ps[:, 0], system.ps[:, 1]
        ax.scatter(x, y, s=sizes, **kwargs)
        ax.set_xlim(-Rmax, Rmax)
        ax.set_ylim(-Rmax, Rmax)
    else:
        x, y, z = system.ps[:, 0], system.ps[:, 1], system.ps[:, 2]
        ax.scatter(x, y, z, s=sizes, **kwargs)
        ax.set_xlim(-Rmax, Rmax)
        ax.set_ylim(-Rmax, Rmax)
        ax.set_zlim(-Rmax, Rmax)
        
    fig, ax = style_plot(ax)
    ax.set_xlabel('X Position (AU)')
    ax.set_ylabel('Y Position (AU)')
    if titl: ax.set_title(titl)
    
    return fig, ax 

def plot_preformance(n_elements, times, labl=None, ax=None, **kwargs):
    """
    Create a properly formatted preformance analysis graph
    
    Parameters
    ----------
    n_elements : array - like
    times : array - like
    labl : string, optional
    ax : matplotlib axis, optional
        If provided, plot on this axis. If None, create new figure
    **kwargs : dict
        Additional plotting parameters. Common options:
        - cmap : str, colormap name (default: 'viridis')
        - s : float, marker size (default: 20)
        - alpha : float, transparency (default: 0.7)
        - figsize : tuple, figure dimensions (only used if ax is None)
        
    Returns
    -------
    fig, ax : matplotlib objects
        Figure and axis used for plotting
        
    """
    
    if ax is None:
        figsize = kwargs.pop('figsize', (10, 8)) 
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    cmap = kwargs.get('cmap', 'coolwarm_r')
    s = kwargs.get('s', 20)
    alpha = kwargs.get('alpha', 0.7)
    
    labl = f'{labl}' if labl else ''
    
    ax.plot(n_elements, times, **kwargs, label=labl)
    fig, ax = style_plot(ax)
    ax.set_xlabel('log Number of Elements')
    ax.set_ylabel('log Time to Complete (s)')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_title('Preformance Analysis')
    if labl: ax.legend()
    
    return fig, ax

def plot_over_time(times, values, labl=None, src=None, titl=None, units=None, eq=None, scale=None, ax=None, **kwargs):
    """
    Create a properly formatted graph of any parameter over time
    
    Parameters
    ----------
    times : array of times
    values : array of values
    labl : labl for legend, optional
    src : string of source names, optional
    titl : title, optional
    units : string of value units, optional
    eq : string of eq for y-axis, optional
    scale : scales y-axes, optional
    ax : matplotlib axis, optional
        If provided, plot on this axis. If None, create new figure
    **kwargs : dict
        Additional plotting parameters. Common options:
        - cmap : str, colormap name (default: 'viridis')
        - s : float, marker size (default: 20)
        - alpha : float, transparency (default: 0.7)
        - figsize : tuple, figure dimensions (only used if ax is None)
        
    Returns
    -------
    fig, ax : matplotlib objects
        Figure and axis used for plotting
        
    """
    
    # Handle ax parameter - create new figure if None
    if ax is None:
        figsize = kwargs.pop('figsize', (8, 6))  # Remove figsize from kwargs
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
    
    # Extract common kwargs with defaults
    cmap = kwargs.get('cmap', 'coolwarm_r')
    s = kwargs.get('s', 20)
    alpha = kwargs.get('alpha', 0.7)
    
    y_labl = f"{f'{scale} ' if scale else ''}{f'{eq} ' if eq else f'{titl} '}{f' ({units})' if units else ''}"
    labl = f'{labl}' if labl else ''
    
    ax.plot(times, values, label=labl, alpha=alpha, **kwargs)    
    fig, ax = style_plot(ax)
    ax.set_xlabel('Time (yrs)')
    ax.set_ylabel(y_labl)
    if scale: ax.set_yscale(scale)
    if titl and src: ax.set_title(f'{src} : {titl}')
    if titl: ax.set_title(titl)
    if labl: ax.legend()
    
    return fig, ax 

def plot_energy_diagnosis(results, src, method, methods=None, scale=None, ax=None, fig=None, **kwargs):
    """
    Creates the energy diagnostic plots
    
    Parameters
    ----------
    results : dataclass
        - Es : Total energies
        - kEs : Kinetic energies
        - pEs : Potential energies
        - v_ratios : virial ratios
        - E_errors : Relative energy errors
        - ts : times
    src : string of source names
    method : string of integration method for labels
    methods : array of methods for legend, optional
    scale : scales y-axis, optional
    ax, fig : matplotlib axis, optional
        If provided, plot on this axis. If None, create new figure
    **kwargs : dict
        Additional plotting parameters. Common options:
        - cmap : str, colormap name (default: 'viridis')
        - s : float, marker size (default: 20)
        - alpha : float, transparency (default: 0.7)
        - figsize : tuple, figure dimensions (only used if ax is None)
        
    Returns
    -------
    fig, ax : matplotlib objects
    Figure and axis used for plotting
    
    """

    if ax is None:
        figsize = kwargs.pop('figsize', (8, 6))
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=figsize)
    else:
        fig = ax[0,0].get_figure() if hasattr(ax, '__getitem__') else ax.get_figure()
    
    N = results.ps.shape[1]
    Eunits = r"$M_\odot\,\mathrm{AU}^2\,\mathrm{yr}^{-2}$"
    v_ratio_EQ = r"$\frac{|2E_K + W|} {|W|}$"
    E_error_EQ = r"$\frac{E_{tot}(t) - E_{tot,0}} {|E_{tot,0}|}$"
    if src is None: src = f'{N}-Body'
    
    comp_colors = {
        'Total Energy':    'C0',
        'Kinetic Energy':  'C1',
        'Potential Energy':'C2',
    }
    ls_map = {'euler':'-', 'RK2':'-.', 'RK4':'--', 'leapfrog':':'}
    
    linestyle = kwargs.pop('linestyle', ls_map.get(method, '-'))
    add_legends = kwargs.pop('add_legends', False)

    fig, ax[0, 1] = plot_over_time(results.ts, results.Es, labl=None, titl='Energy Components', src=src,
                                   units=Eunits, scale=None, ax=ax[0, 1], color=comp_colors['Total Energy'], linestyle=linestyle, **kwargs)
    fig, ax[0, 1] = plot_over_time(results.ts, results.kEs, labl=None, titl='Energy Components', src=src,
                                   units=Eunits, scale=None, ax=ax[0, 1], color=comp_colors['Kinetic Energy'], linestyle=linestyle, **kwargs)
    fig, ax[0, 1] = plot_over_time(results.ts, results.pEs, labl=None, titl='Energy Components', src=src,
                                   units=Eunits, scale=None, ax=ax[0, 1], color=comp_colors['Potential Energy'], linestyle=linestyle, **kwargs)
    
    fig, ax[1, 0] = plot_over_time(results.ts, results.E_errors, titl='Relative Energy Errors', labl=method, eq=E_error_EQ, src=src, 
                                   scale=scale, ax=ax[1, 0], linestyle=linestyle, **kwargs)
    fig, ax[1, 1] = plot_over_time(results.ts, results.v_ratios, titl='Virial Ratios', labl=method, eq=v_ratio_EQ, src=src, 
                                   scale=scale, ax=ax[1, 1], linestyle=linestyle, **kwargs)

    if add_legends: style_energy_comps(ax[0, 1], methods)
    
    return fig, ax

def style_energy_comps(ax, used_methods=None):
    """
    Add standardized legends to the Energy Components subplot and only include
    integrators that were actually plotted (based on linestyles present on `ax`).

    Parameters
    ----------
    ax : matplotlib axis
        Axis that contains the energy component lines
    used methods : array of strings of methods used
    
    """

    comp_colors = {
        'Total Energy':     'C0',
        'Kinetic Energy':   'C1',
        'Potential Energy': 'C2',
    }
    
    ls_map = {'euler':'-', 'RK2':'-.', 'RK4':'--', 'leapfrog':':'}

    comp_handles = [
        Line2D([0],[0], color=comp_colors['Total Energy'], lw=2, label='Total'),
        Line2D([0],[0], color=comp_colors['Kinetic Energy'], lw=2, label='Kinetic'),
        Line2D([0],[0], color=comp_colors['Potential Energy'], lw=2, label='Potential'),
    ]
    
    comp_legend = ax.legend(handles=comp_handles, title='Energy component', loc='upper right')

    if used_methods is None:
        present_styles = {line.get_linestyle() for line in ax.get_lines()}
        used_methods = [name for name, style in ls_map.items() if style in present_styles]

    if used_methods:
        style_handles = [
            Line2D([0],[0], color='k', lw=2, linestyle=ls_map[name], label=name)
            for name in used_methods if name in ls_map
        ]
        ax.add_artist(comp_legend)
        ax.legend(handles=style_handles, title='Integrator', loc='center right')

    return ax

def histogram(values, xlabl=None, titl=None, ax=None, kroupa=None, **kwargs):
    """
    Create a properly formatted histogram
    
    Parameters
    ----------
    values : array of values
    xlabl : xlabel string
    titl : title string
    ax : matplotlib axis, optional
        If provided, plot on this axis. If None, create new figure
    kroupa : Boolean, set True if you want to overplot kroupa lines
    **kwargs : dict
        Additional plotting parameters. Common options:
        - cmap : str, colormap name (default: 'viridis')
        - s : float, marker size (default: 20)
        - alpha : float, transparency (default: 0.7)
        - figsize : tuple, figure dimensions (only used if ax is None)
        
    Returns
    -------
    fig, ax : matplotlib objects
        Figure and axis used for plotting
        
    """
    
    if ax is None:
        figsize = kwargs.pop('figsize', (8, 6))  
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
    
    cmap = kwargs.get('cmap', 'coolwarm_r')
    s = kwargs.get('s', 20)
    alpha = kwargs.get('alpha', 0.7)
    
    xlabl = f'{xlabl}' if xlabl else f''
    
    vals = np.asarray(values).ravel()
    vals = vals[np.isfinite(vals) & (vals > 0)]
    bins = np.geomspace(vals.min(), vals.max(), 50)
    
    ax.hist(values, bins=bins, density=True, **kwargs, linewidth=1, edgecolor='k')    
    fig, ax = style_plot(ax)
    ax.set_xscale('log')              
    ax.set_yscale('log')
    ax.set_xlabel(xlabl)
    ax.set_ylabel('Number Density')
    if titl: ax.set_title(titl)
    
    if kroupa:
        
        ax.axvline(0.5, ls='--', c='red', lw=3, label=fr'0.5 M$_{{\odot}}$')

        dens, edges = np.histogram(values, bins=bins, density=True)
        idx = np.searchsorted(edges, 0.5, side='right') - 1
        idx = np.clip(idx, 0, len(dens) - 1)
        y0 = dens[idx] if np.isfinite(dens[idx]) and dens[idx] > 0 else np.nanmax(dens)

        x1 = np.geomspace(edges[0], 0.5, 200)
        y1 = y0 * (x1 / 0.5) ** (-1.3)

        x2 = np.geomspace(0.5, edges[-1], 200)
        y2 = y0 * (x2 / 0.5) ** (-2.3)

        ax.plot(x1, y1, lw=1.5, label=fr'$m^-{1.3}$')
        ax.plot(x2, y2, lw=1.5, label=fr'$m^-{2.3}$')
        
        ax.legend()
    
    return fig, ax 

def plot_density_profile(radials, a, labl=None, titl=None, ax=None, **kwargs):
    """
    Plot radial mass-density and a Plummer overlay (scaled to match at r=a).

    Parameters
    ----------
    radials : (radii, densities) 
    a : Plummer scale radius (AU)
    labl : label for legend, optional
    titl : label for title, optinal
    ax : matplotlib axis, optional
        If provided, plot on this axis. If None, create new figure
    **kwargs : dict
        Additional plotting parameters. Common options:
        - cmap : str, colormap name (default: 'viridis')
        - s : float, marker size (default: 20)
        - alpha : float, transparency (default: 0.7)
        - figsize : tuple, figure dimensions (only used if ax is None)
        
    Returns
    -------
    fig, ax : matplotlib objects
        Figure and axis used for plotting
    """
    rs, ds = radials
    rs = np.asarray(rs).ravel()
    ds = np.asarray(ds).ravel()

    if ax is None:
        figsize = kwargs.pop('figsize', (7, 5))
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    rplot = np.geomspace(rs.min(), rs.max(), 400)
    model = (1.0 + (rplot / a)**2)**(-2.5)

    r_ref = np.clip(a, rs.min(), rs.max())
    y_ref_data  = np.interp(r_ref, rs, ds)
    y_ref_model = (1.0 + (r_ref / a)**2)**(-2.5)
    scale = y_ref_data / y_ref_model if y_ref_model > 0 else 1.0
    
    ax.plot(rs, ds, '.', ms=3, lw=1.2, label=labl, **kwargs)
    ax.plot(rplot, scale * model, '--', lw=1.5, label=r'Plummer: $\rho\propto(1+r^2/a^2)^{-5/2}$')
    fig, ax = style_plot(ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('r (AU)')
    ax.set_ylabel(r'$\rho\ \left(M_\odot\,\mathrm{AU}^{-3}\right)$')
    if titl: ax.set_title(titl)
    ax.legend()
    
    return fig, ax

def many_snapshots(results_x, a, ms=None, ax=None, fig=None, **kwargs):
    """

    Parameters
    ----------
    results_x : Results dataclass
    a : half mass radius (AU)
    ms : array of masses (Msun), optional
    ax : matplotlib axis, optional
    fig : matplotlib figure, optional
    **kwargs : dict
        Additional plotting parameters. Common options:
        - cmap : str, colormap name (default: 'viridis')
        - s : float, marker size (default: 20)
        - alpha : float, transparency (default: 0.7)
        - figsize : tuple, figure dimensions (only used if ax is None)
    
    Returns
    -------
    fig, ax : matplotlib objects
        Figure and axis used for plotting
    """
    snaps, N, d = results_x.ps.shape

    if snaps % 2 == 0:
        rows = 2
        cols = snaps // 2
    elif snaps % 3 == 0:
        rows = 3
        cols = snaps // 3
    else:
        rows = 2
        cols = (snaps + 1) // 2

    if ax is None:
        if d == 3:
            fig = plt.figure(figsize=(cols*4, rows*4))
            ax = np.empty((rows, cols), dtype=object)
            for k in range(rows * cols):
                i, j = divmod(k, cols)
                ax[i, j] = fig.add_subplot(rows, cols, k + 1, projection='3d')
        else:
            fig, ax = plt.subplots(nrows=rows, ncols=cols, figsize=(rows*4, cols*4))
    else:
        fig = ax.get_figure()

    Rmax = 2*a

    class _S: pass

    for s in range(snaps):
        i, j = divmod(s, cols)
        sys_s = _S(); sys_s.ms = np.ones((N, 1)); sys_s.ps = results_x.ps[s]; sys_s.dimensions = d

        frac = s / (snaps - 1) if snaps > 1 else 0.0
        titl = rf"$t/t_{{\rm sim}}$ = {frac:.2f}"

        fig, ax[i, j] = plot_snapshot(sys_s, a, ms, titl=titl, ax=ax[i, j], **kwargs)

        if d == 3:
            ax[i, j].set_xlim(-Rmax, Rmax); ax[i, j].set_ylim(-Rmax, Rmax); ax[i, j].set_zlim(-Rmax, Rmax)
        else:
            ax[i, j].set_aspect('equal', 'box')
            ax[i, j].set_xlim(-Rmax, Rmax); ax[i, j].set_ylim(-Rmax, Rmax)

    for k in range(snaps, rows * cols):
        i, j = divmod(k, cols)
        ax[i, j].axis('off')

    return fig, ax

def make_orbit_gif(results, out='orbits.gif', titl=None):
    """
    Parameters
    ----------
    results : Results dataclass
    out : filename, optional
    titl : title string, optional
    **kwargs : dict
        Additional plotting parameters. Common options:
        - cmap : str, colormap name (default: 'viridis')
        - s : float, marker size (default: 20)
        - alpha : float, transparency (default: 0.7)
        - figsize : tuple, figure dimensions (only used if ax is None)
    
    Returns
    -------
    nothing, saves plot
    """
    
    plane = 'xy'
    s = 14
    color = 'k'
    fps = 20
    dpi = 100
    
    i, j = {'xy': (0,1), 'xz': (0,2), 'yz': (1,2)}.get(plane, (0,1))
    S, N, d = results.ps.shape

    fig, ax = plt.subplots(figsize=(6,6))
    plot_orbit(results, labl='body', titl=titl, ax=ax, alpha=0.35, linewidth=0.8)

    sc = ax.scatter(results.ps[0,:,i], results.ps[0,:,j], s=s, c=color, alpha=0.9)
    txt = ax.text(0.02, 0.98, '', transform=ax.transAxes, va='top')

    def update(k):
        sc.set_offsets(np.column_stack((results.ps[k,:,i], results.ps[k,:,j])))
        frac = k / (S-1) if S>1 else 0.0
        txt.set_text(rf'$t/t_{{\rm sim}}={frac:.2f}$')
        return sc, txt

    anim = FuncAnimation(fig, update, frames=S, interval=1000/fps, blit=False)
    anim.save('gifs/'+out, writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
