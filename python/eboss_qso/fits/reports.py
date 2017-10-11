from pyRSD.rsdfit import FittingDriver
from pyRSD.rsdfit.results import LBFGSResults, EmceeResults
from pyRSD.rsdfit.analysis import BestfitParameterSet, tex_names
from pyRSD.rsdfit.plotly import *
from glob import glob
import os
from jinja2 import Environment, FileSystemLoader

from eboss_qso.measurements import get_hashkeys
from eboss_qso import EBOSS_FITS
from .driver import mkdir_p

def _summarize_fit(info,  d):
    """
    Internal function to summarize the fit parameters by yielding HTML list
    items with each parameter.

    Parameters
    ----------
    info : dict
        this is the dictionary of hash keys for the fit
    d : FittingDriver
        the fitting driver holding the fit results

    Yields
    ------
    str :
        HTML list group item holding parameter key/values
    """
    sample = "NGC" if 'ngc' in d.data.statistics[0] else "SGC"

    out = []
    out.append(("sample", sample))
    out.append((r"$P_0^\mathrm{FKP}$", "%.2e" %info['P0_FKP']))
    out.append((r"$(z_\mathrm{min}, \ z_\mathrm{max})$", (info['zmin'], info['zmax'])))
    out.append(("redshift weighted", info['z-weighted']))
    out.append(("$p$", info['p']))
    out.append(("statistics", info['stats']))
    out.append(("free parameters", info['vary']))
    out.append(("$k$ range", d.data.fitting_range))
    out.append(("burn-in steps", d.results.burnin))

    for args in out:
        yield r"<li class='list-group-item'>%s = %s</li>" % args

def _load_fit(dirname, burnin=None):
    """
    Internal function to load a QSO fit from the specified directory.

    The driver is initialized from the most recent results file in the
    input directory.

    Returns
    -------
    FittingDriver :
        the driver object holding the result
    """
    pattern = os.path.join(dirname, '*npz')
    files = sorted(glob(pattern), key=os.path.getmtime, reverse=True)
    if not len(files):
        raise ValueError("no result files found in directory '%s'" %dirname)
    d = FittingDriver.from_directory(dirname, results_file=files[0])

    if burnin is None:
        burnin = int(0.5*d.results.iterations)
    d.results.burnin = burnin
    d.set_fit_results()
    return d

def generate_redshift_summary(dirname, output, burnin=None):
    """
    For a given fit result directory, generate an HTML report.

    Parameters
    ----------
    dirname : str
        the name of the directory holding the fit results
    output : str
        the output file name
    burnin : int, optional
        the burn-in period
    """
    # get the template path/dir
    template_dir = os.path.join(EBOSS_FITS, "html-templates")
    template_file = os.path.join(template_dir, "report-template.html")

    # load the driver and bestfit results
    d = _load_fit(dirname, burnin=burnin)
    if isinstance(d.results, EmceeResults):
        bestfit = BestfitParameterSet.from_mcmc(d.results)
    else:
        bestfit = BestfitParameterSet.from_nlopt(d.results)
    bestfit['tex_name'] = [tex_names.get(name, name) for name in bestfit.index]

    # get the hashkeys
    info = get_hashkeys(dirname, None)

    # make the fit summary
    fit_summary = "<ul class='list-group'>"
    for s in _summarize_fit(info, d):
        fit_summary += s + '\n'
    fit_summary += "</ul>"

    # get the table
    table = bestfit.to_ipynb().to_html()
    table = table.replace('<table border="1" class="dataframe">','<table class="table table-striped">')

    # plot theory/data comparison
    labels = []
    for stat in info['stats']:
        if stat == "P0_sysfree":
            labels.append(r"$P_0 + 2/5 P_2$")
        else:
            labels.append(r"$%s_%s$" %(stat[0], stat[1]))

    fig = plot_fit_comparison(d, labels=labels)
    div1 = py.plot(fig, output_type='div', include_plotlyjs=False)

    # plot the triangle
    div2 = ""
    if isinstance(d.results, EmceeResults):
        N = len(d.results.free_names)
        width = 250*N; height = 250*N
        fig = plot_triangle(d.results, params=d.results.free_names, thin=5, width=width, height=height)
        div2 = py.plot(fig, output_type='div', include_plotlyjs=False, image_width=width, image_height=height)

    # plot the traces
    fig = plot_traces(d.results, *d.results.free_names, burnin=0, max_walkers=20)
    div3 = py.plot(fig, output_type='div', include_plotlyjs=False)

    # render and save!
    jinja_env = Environment(loader=FileSystemLoader(template_dir))
    tpl = jinja_env.get_template(os.path.basename(template_file))
    html_file = tpl.render(fit_summary=fit_summary, summary_table=table, summary_plot=div1, histogram=div2, traces=div3)

    with open(output, 'w') as ff:
        ff.write(html_file)


def generate_fit_report(dirname, output, burnin=None):
    """
    For a given fit result directory, generate an HTML report.

    Parameters
    ----------
    dirname : str
        the name of the directory holding the fit results
    output : str
        the output file name
    burnin : int, optional
        the burn-in period
    """
    # get the template path/dir
    template_dir = os.path.join(EBOSS_FITS, "html-templates")
    template_file = os.path.join(template_dir, "report-template.html")

    # load the driver and bestfit results
    d = _load_fit(dirname, burnin=burnin)
    if isinstance(d.results, EmceeResults):
        bestfit = BestfitParameterSet.from_mcmc(d.results)
    else:
        bestfit = BestfitParameterSet.from_nlopt(d.results)
    bestfit['tex_name'] = [tex_names.get(name, name) for name in bestfit.index]

    # get the hashkeys
    info = get_hashkeys(dirname, None)

    # make the fit summary
    fit_summary = "<ul class='list-group'>"
    for s in _summarize_fit(info, d):
        fit_summary += s + '\n'
    fit_summary += "</ul>"

    # get the table
    table = bestfit.to_ipynb().to_html()
    table = table.replace('<table border="1" class="dataframe">','<table class="table table-striped">')

    # plot theory/data comparison
    labels = []
    for stat in info['stats']:
        if stat == "P0_sysfree":
            labels.append(r"$P_0 + 2/5 P_2$")
        else:
            labels.append(r"$%s_%s$" %(stat[0], stat[1]))

    fig = plot_fit_comparison(d, labels=labels)
    div1 = py.plot(fig, output_type='div', include_plotlyjs=False)

    # plot the triangle
    div2 = ""
    if isinstance(d.results, EmceeResults):
        N = len(d.results.free_names)
        width = 250*N; height = 250*N
        fig = plot_triangle(d.results, params=d.results.free_names, thin=5, width=width, height=height)
        div2 = py.plot(fig, output_type='div', include_plotlyjs=False, image_width=width, image_height=height)

    # plot the traces
    fig = plot_traces(d.results, *d.results.free_names, burnin=0, max_walkers=20)
    div3 = py.plot(fig, output_type='div', include_plotlyjs=False)

    # render and save!
    jinja_env = Environment(loader=FileSystemLoader(template_dir))
    tpl = jinja_env.get_template(os.path.basename(template_file))
    html_file = tpl.render(fit_summary=fit_summary, summary_table=table, summary_plot=div1, histogram=div2, traces=div3)

    with open(output, 'w') as ff:
        ff.write(html_file)

def _get_params_from_model_string(model):
    """
    Internal function to return a string of Latex'ed parameter names
    based on the input model string.
    """
    from pyRSD.rsdfit.analysis import tex_names
    out = []
    model = model.split('-')
    if 'basemodel' in model:
        out += [tex_names[p] for p in ['b1', 'sigma_fog']]
    if 'fs8' in model:
        out += [tex_names[p] for p in ['f', 'sigma8_z']]
    if 'alphas' in model:
        out += [tex_names[p] for p in ['alpha_par', 'alpha_perp']]
    if 'fnl' in model:
        out += [tex_names['f_nl']]

    return "[" + ", ".join(out) + "]"

def generate_toc(kind):
    """
    Generate the table of contents for the reports, using the structure of
    the HTML reports directory tree.
    """
    # either REPORT_DIR/data or REPORT_DIR/mocks/ezmock
    if kind == 'data':
        subpath = kind
    else:
        subpath = os.path.join('mocks', kind)

    # get the template path/dir
    template_dir = os.path.join(EBOSS_FITS, "html-templates")
    template_file = os.path.join(template_dir, "index-template.html")

    # where we start looking for reports
    home_dir = os.path.join(EBOSS_FITS, 'reports', subpath)
    versions = sorted(glob(os.path.join(home_dir, '*')), reverse=True)

    out = ""
    for v in versions:
        v = os.path.split(v)[-1]
        if v == 'index.html': continue
        out += f"<h1 id='%s'>Version %s</h1>\n" %(v[1:], v[1:])

        d0 = os.path.join(home_dir, v)
        kranges = sorted(glob(os.path.join(d0, '*')))
        for krange in kranges:
            krange = os.path.split(krange)[-1]
            kmin,kmax = krange.split('-')
            out += r"<h2 id='%s'>Fitting Range: (%s, %s)</h2>" %(krange, kmin, kmax)
            out += "\n"

            d1 = os.path.join(d0, krange)
            models = sorted(glob(os.path.join(d1, '*')))
            for model in models:
                model = os.path.split(model)[-1]
                params = _get_params_from_model_string(model)
                out += r"<h3 id='%s'>Parameter Set: %s</h3>" % (model, params)
                out += "\n"

                d2 = os.path.join(d1, model)
                zranges = sorted(glob(os.path.join(d2, '*')))
                for zrange in zranges:
                    zrange = os.path.split(zrange)[-1]
                    zmin,zmax = zrange.split('-')
                    out += r"<h4 id='%s'>Redshift Range: (%s, %s)</h2>" %(zrange, zmin, zmax)
                    out += "\n"

                    d3 = os.path.join(d2, zrange)
                    fits = sorted(glob(os.path.join(d3, '*')))
                    out += "<ul>\n"
                    for fit in fits:
                        relpath = os.path.join(os.path.relpath(fit, home_dir), 'report.html')
                        fit = os.path.split(fit)[-1]
                        out += "<li><a href='%s'>%s</a></li>\n" %(relpath, fit)
                    out += "</ul>\n"

    # render and save!
    jinja_env = Environment(loader=FileSystemLoader(template_dir))
    tpl = jinja_env.get_template(os.path.basename(template_file))
    html_file = tpl.render(index=out)

    output = os.path.join(home_dir, 'index.html')
    print('saving new TOC: %s' %output)
    with open(output, 'w') as ff:
        ff.write(html_file)


def generate_all(kind, dirpaths=[], overwrite=False):
    """
    Generate the HTML reports, optionally overwriting existing reports.
    """
    if kind == 'data':
        subpath = kind
    else:
        subpath = os.path.join('mocks', kind)
    results_dir = os.path.join(EBOSS_FITS, 'results', subpath)
    reports_dir = os.path.join(EBOSS_FITS, 'reports', subpath)

    def _generate(dirpath):

        # find the results file
        pattern = os.path.join(dirpath, '*npz')
        result_file = sorted(glob(pattern), key=os.path.getmtime, reverse=True)[0]
        mtime = os.path.getmtime(result_file)

        # get the matching report directory
        relpath = os.path.relpath(dirpath, results_dir)
        this_report_dir = os.path.join(reports_dir, relpath)
        r = os.path.join(this_report_dir, 'report.html')
        if os.path.isdir(this_report_dir) and os.path.isfile(r):
            if os.path.getmtime(r) >= mtime and not overwrite:
                return 0

        # need to make the report
        if not os.path.isdir(this_report_dir):
            mkdir_p(this_report_dir)
        print("generating report for %s..." % relpath)
        generate_fit_report(dirpath, r)
        return 1

    # do specific reports
    N = 0
    if len(dirpaths):
        for dirpath in dirpaths:
            N += _generate(os.path.abspath(dirpath))
    # do all out-of-date reports
    else:
        # walk the full directory path
        for dirpath, dirnames, filenames in os.walk(results_dir):

            # this is a fit result directory
            if all(f in filenames for f in ['params.dat', 'hashinfo.json']):
                if any(f.endswith('.npz') for f in filenames):
                    N += _generate(dirpath)

    return N

def _generate_all():
    """
    The console script to generate HTML scripts, which reads in command-line
    arguments from the user.
    """
    import argparse

    desc = 'generate HTML reports for eBOSS QSO fits'
    parser = argparse.ArgumentParser(description=desc)

    h = 'the kind of fits, i.e., data, ezmock, etc'
    parser.add_argument('kind', choices=['data', 'ezmock'], type=str, help=h)

    h = 'directories to generate reports for'
    parser.add_argument('--dirpaths', nargs='*', default=[], type=str, help=h)

    h = 'whether to overwrite any existing reports'
    parser.add_argument('--overwrite', action='store_true', help=h)

    ns = parser.parse_args()

    # generate all reports
    updated = generate_all(ns.kind, dirpaths=ns.dirpaths, overwrite=ns.overwrite)

    # generate the new TOC
    if updated:
        generate_toc(ns.kind)