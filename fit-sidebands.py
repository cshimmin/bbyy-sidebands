#!/usr/bin/env python

import ROOT as r
from glob import glob
import numpy as np
import sys
import os
import array as arr

# myybb range to fit for the lowmass analysis
MASS_RANGE_LOW = (245, 610)
NBINS_LOW = 56

# myybb range to fit for the highmass analysis
MASS_RANGE_HIGH = (335, 1140)
NBINS_HIGH = 50

BTAG_CAT_DEFAULT = 0

LABEL_POS = (0.6,0.7,0.8,0.9)


# Define a bunch of different ways to slice the myy bins.
# You can select which one to use w/ the command line option
BOUNDARIES_DEFAULT = 'sr-window-1' 
BOUNDARIES = {
        # sr-window-1 is defined such that exactly one slice
        # fits in the signal region, sr-window-2 has two slices
        # in the SR, etc.
        'sr-window-1': (
            [
                (111, 120.4),
                (129.8, 139.2),
                (139.2, 148.6),
                (148.6, 158),
            ],
            [
                (120.4, 129.8)
            ]),
        'sr-window-2': (
            [
                (106.3, 111.0),
                (111.0, 115.7),
                (115.7, 120.4),
                (129.8, 134.5),
                (134.5, 139.2),
                (139.2, 143.9),
                (143.9, 148.6),
                (148.6, 153.3),
                (153.3, 158.0),
            ],
            [
                (120.4, 125.1),
                (125.1, 129.8),
            ]),
        'sr-window-4': (
            [
                (106.30, 108.65),
                (108.65, 111.00),
                (111.00, 113.35),
                (113.35, 115.70),
                (115.70, 118.05),
                (118.05, 120.40),
                (129.80, 132.15),
                (132.15, 134.50),
                (134.50, 136.85),
                (136.85, 139.20),
                (139.20, 141.55),
                (141.55, 143.90),
                (143.90, 146.25),
                (146.25, 148.60),
                (148.60, 150.95),
                (150.95, 153.30),
                (153.30, 155.65),
                (155.65, 158.00),
            ],
            [
                (120.40, 122.75),
                (122.75, 125.10),
                (125.10, 127.45),
                (127.45, 129.80),
                #(121,123),
                #(123,125),
                #(125,127),
                #(127,129),
            ]),
            'sb-low-1': (
            [
                (105.0,120.4),
                (129.8,145.2),
            ],
            []),
            'sb-low-2': (
            [
                (105.0,112.7),
                (112.7,120.4),
                (129.8,137.5),
                (137.5,145.2),
                (145.2,152.9),
            ],
            []),
        }

def up(interact=True):
    r.gROOT.FindObject("c1").Update()
    if interact:
        raw_input("Press enter")

def read_tree(t):
    data = []
    for evt in t:
        myy = getattr(t, "HGamEventInfoAuxDyn.m_yy")
        myybb = getattr(t, "HGamEventInfoAuxDyn.yybb_lowMass_m_yybb")
        #cat = getattr(t, "HGamEventInfoAuxDyn.yybb_bTagCat")

        #data.append((myy/1e3,myybb/1e3,cat))
        data.append((myy/1e3,myybb/1e3))
    #return np.array(data, dtype=[('myy', float), ('myybb', float), ('cat', int)])
    return np.array(data, dtype=[('myy', float), ('myybb', float)])

MODEL_OPTIONS = ('novo', 'landau', 'mlandau', 'gamma', 'exp', 'epoly2', 'chase-1')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--binning", choices=sorted(BOUNDARIES.keys()), default=BOUNDARIES_DEFAULT, help="myy binning scheme to use (default: %s)"%BOUNDARIES_DEFAULT)
    parser.add_argument("--btag", type=int, choices=(0,1,2), default=BTAG_CAT_DEFAULT, help="the btag category to use.")
    parser.add_argument("--highmass", action="store_true", help="use high mass selection")
    parser.add_argument("--data-path", default="data-skim", help="the path containing the (skimmed) MxAOD data (default: ./data-skim)")
    parser.add_argument("--sr", action="store_true", help="use the 0tag signal region")
    parser.add_argument("--out", default=".", help="directory in which to output images and parameters")
    parser.add_argument("--fastforward", action="store_true", help="don't pause after individual sideband slice fits")
    parser.add_argument("--model", default=MODEL_OPTIONS[0], choices=MODEL_OPTIONS, help="BG model to use. (default: %s)"%MODEL_OPTIONS[0])
    parser.add_argument("--tail", type=float, help="set the tail param to a specific value")
    parser.add_argument("--no-interact", action="store_true", help="skip interactive prompts.")
    args = parser.parse_args()

    if args.sr and not args.btag == 0:
        print >> sys.stderr, "ERROR: Can only do SR fits in 0-tag data."
        sys.exit(1)

    try:
        os.makedirs(args.out)
    except OSError as e:
        import errno
        if e.errno == errno.EEXIST:
            pass
        else:
            raise

    if args.highmass:
        MASS_RANGE = MASS_RANGE_HIGH
        NBINS = NBINS_HIGH
    else:
        MASS_RANGE = MASS_RANGE_LOW
        NBINS = NBINS_LOW

    ds_boundaries, ds_boundaries_sr = BOUNDARIES[args.binning]

    t0 = r.TChain("CollectionTree")
    map(t0.Add, glob(os.path.join(args.data_path, "*.root")))

    if t0.GetEntries() == 0:
        print >> sys.stderr, "ERROR: Found 0 events in path:", args.data_path
        sys.exit(1)

    masscat = "highMass" if args.highmass else "lowMass"
    selection = "HGamEventInfoAuxDyn.yybb_%s_cutFlow==4 && HGamEventInfoAuxDyn.yybb_bTagCat==%d"%(masscat, args.btag)

    t = t0.CopyTree(selection)

    data = read_tree(t)

    #data = data[data['cat']==args.btag]

    print "got %d events in 0tag" % len(data)

    w = r.RooWorkspace("w")

    myybb = w.factory("myybb[%g,%g]"%MASS_RANGE)
    myybb.setUnit("GeV")
    myybb.setBins(NBINS)
    myybb_ = r.RooArgSet(myybb)

    w.factory("win_low[105]")
    w.factory("win_hi[160]")

    ds_lo = r.RooDataSet("ds_lo","ds_lo",myybb_)
    ds_hi = r.RooDataSet("ds_hi","ds_hi",myybb_)

    datasets_sr = []
    if args.sr:
        #t_sr_ = r.TChain("CollectionTree")
        #map(t_sr_.Add, glob("data-0tagSR/*.root"))

        selection_sr = "HGamEventInfoAuxDyn.yybb_%s_cutFlow>4 && HGamEventInfoAuxDyn.yybb_bTagCat==0"%masscat
        t_sr = t0.CopyTree(selection_sr)

        data_sr = read_tree(t_sr)

        for i in xrange(len(ds_boundaries_sr)):
            datasets_sr.append(r.RooDataSet("ds_sr%d"%i,"ds_sr%d"%i,myybb_))

        for xyy, xyybb in data_sr:
            if xyybb < MASS_RANGE[0] or xyybb > MASS_RANGE[1]: continue
            myybb.setVal(xyybb)
            
            for i,(xlo,xhi) in enumerate(ds_boundaries_sr):
                if xyy>=xlo and xyy<xhi:
                    datasets_sr[i].add(myybb_)

    datasets = []
    boundary_averages = [0]*len(ds_boundaries)
    for i in xrange(len(ds_boundaries)):
        datasets.append(r.RooDataSet("ds%d"%i,"ds%d"%i,myybb_))

    for xyy, xyybb in data:
        if xyybb < MASS_RANGE[0] or xyybb > MASS_RANGE[1]: continue
        myybb.setVal(xyybb)
        if xyy < 125:
            ds_lo.add(myybb_)
        else:
            ds_hi.add(myybb_)

        for i,(xlo,xhi) in enumerate(ds_boundaries):
            if xyy>=xlo and xyy<xhi:
                datasets[i].add(myybb_)
                boundary_averages[i] += xyy

    for i in xrange(len(ds_boundaries)):
        boundary_averages[i] = boundary_averages[i] / datasets[i].sumEntries()

    if args.model == 'novo':
        pdf = w.factory("Novosibirsk::bg_pdf(myybb, bg_peak[200,400], bg_width[0,200], bg_tail[-2.5,-0.01])")
        ndof = 3
        if not args.tail is None:
            w.obj("bg_tail").setVal(args.tail)
            w.obj("bg_tail").setConstant(True)
            ndof = 2
        param_names = ['bg_peak','bg_width','bg_tail']
    elif args.model == 'landau':
        pdf = w.factory("Landau::bg_pdf(myybb, bg_peak[200,350], bg_width[0,200])")
        ndof = 2
        param_names = ['bg_peak','bg_width']
    elif args.model == 'mlandau':
        pdf = w.factory("Landau::bg_pdf(myybb, bg_peak[200,350], expr::bg_scale('bg_p0 + bg_p1*(myybb-250)', bg_p0[0,200], bg_p1[-10,10], myybb))")
        ndof = 3
        param_names = ['bg_peak', 'bg_p0', 'bg_p1']
    elif args.model == 'gamma':
        pdf = w.factory("Gamma::bg_pdf(myybb, bg_gamma[0,100], bg_beta[0,100], bg_mu[0,300])")
        ndof = 3
        param_names = ['bg_mu', 'bg_gamma', 'bg_beta']
    elif args.model == 'mgamma':
        w.factory("expr::bg_gamma('bg_y0 + bg_y1*(myybb-200)', bg_y0[0,10], bg_y1[-1,1])")
        w.factory("expr::bg_beta('bg_b0 + bg_b1*(myybb-200)', bg_b0[0,100], bg_b1[-1,1])")
        pdf = w.factory("Gamma::bg_pdf(myybb, bg_gamma, bg_beta, bg_mu[0,300])")
        ndof = 5
        param_names = ['bg_mu', 'bg_y0', 'bg_y1', 'bg_b0', 'bg_b1']
    elif args.model == 'exp':
        pdf = w.factory("Exponential::bg_pdf(myybb, bg_slope[-0.1,0])")
        ndof = 1
        param_names = ['bg_slope']
    elif args.model == 'epoly2':
        pdf = w.factory("EXPR::bg_pdf('exp(bg_p0*myybb*(1 + bg_p1*myybb/100.))', myybb, bg_p0[-0.1,-1e-5], bg_p1[-0.5,0.5])")
        ndof = 2
        param_names = ['bg_p0', 'bg_p1']
    elif args.model == 'chase-1':
        pdf = w.factory("EXPR::bg_pdf('(exp(-a*win_low - myybb*(b0 + b1*win_low/100.)) - exp(-a*win_hi - myybb*(b0 + b1*win_hi/100.)))/(a + b1*myybb/100.)', myybb, a[-1,0], b0[0,1], b1[-10,10], win_low, win_hi)")
        #pdf = w.factory("EXPR::bg_pdf('exp(-a*win_low - myybb*(b0 + b1*win_low)) ', myybb, a[-1,0], b0[-1,0], b1[-10,10], win_low, win_hi)")
        #param_names = ['a','b0','b1']
        #w.obj("a").setVal(-2.59e-2)
        #w.obj("a").setConstant(True)
        ndof = 3
        param_names = ['a','b0','b1']

    w.defineSet("params",",".join(param_names))
    params = w.set("params")

    w.saveSnapshot("init", params)

    frame = myybb.frame()
    ds_lo.plotOn(frame)

    w.loadSnapshot("init")
    pdf.fitTo(ds_lo)
    pdf.plotOn(frame)
    chi2 = frame.chiSquare(ndof)
    frame.SetTitle("%d-tag LOW sideband [%d events]"%(args.btag, ds_lo.sumEntries()))
    frame.Draw()
    txt = r.TPaveLabel(*(LABEL_POS+("#chi^{2}/ndof = %0.2f"%chi2,"NDC")))
    txt.SetFillStyle(0)
    txt.SetBorderSize(0)
    txt.SetTextSize(0.25)
    txt.Draw()
    up((not args.no_interact) and (not args.fastforward))
    r.gROOT.FindObject("c1").SaveAs(os.path.join(args.out,"%dtag_low_sideband.pdf"%args.btag))

    frame = myybb.frame()
    ds_hi.plotOn(frame)

    w.loadSnapshot("init")
    pdf.fitTo(ds_hi)
    pdf.plotOn(frame, r.RooFit.LineColor(r.kRed))
    chi2 = frame.chiSquare(ndof)
    frame.SetTitle("%d-tag HIGH sideband [%d events]"%(args.btag, ds_hi.sumEntries()))
    frame.Draw()
    txt = r.TPaveLabel(*(LABEL_POS+("#chi^{2}/ndof = %0.2f"%chi2,"NDC")))
    txt.SetFillStyle(0)
    txt.SetBorderSize(0)
    txt.SetTextSize(0.25)
    txt.Draw()
    up((not args.no_interact) and (not args.fastforward))
    r.gROOT.FindObject("c1").SaveAs(os.path.join(args.out,"%dtag_high_sideband.pdf"%args.btag))

    colors = [r.kRed, r.kOrange+1, r.kGreen, r.kBlue, r.kViolet]


    fit_vals = []
    fit_errs = []
    fit_vals_sr = []
    fit_errs_sr = []

    x_vals = []
    x_errs = []
    x_vals_sr = []
    x_errs_sr = []

    if args.sr:
        ds_boundaries_all = ds_boundaries + ds_boundaries_sr
        datasets_all = datasets + datasets_sr
    else:
        ds_boundaries_all = ds_boundaries
        datasets_all = datasets

    frame = myybb.frame()
    frame_sr = myybb.frame()
    for i,(xlo,xhi) in enumerate(ds_boundaries_all):
        is_sr = datasets_all[i] in datasets_sr

        x0 = (0.5*(xhi+xlo))
        dx = (0.5*(xhi-xlo))

        w.obj("win_low").setVal(xlo)
        w.obj("win_hi").setVal(xhi)

        w.loadSnapshot("init")
        result = pdf.fitTo(datasets_all[i], r.RooFit.Save(True))

        fit_res = [w.obj(p).getVal()   for p in param_names]
        fit_err = [w.obj(p).getError() for p in param_names]

        if is_sr:
            pdf.plotOn(frame_sr, r.RooFit.LineColor(colors[i%len(colors)]))
            x_vals_sr.append(x0)
            x_errs_sr.append(dx)
            fit_vals_sr.append(fit_res)
            fit_errs_sr.append(fit_err)
        else:
            pdf.plotOn(frame, r.RooFit.LineColor(colors[i%len(colors)])) 
            x_vals.append(x0)
            x_errs.append(dx)
            fit_vals.append(fit_res)
            fit_errs.append(fit_err)

        frame_tmp = myybb.frame()
        datasets_all[i].plotOn(frame_tmp)
        pdf.plotOn(frame_tmp, r.RooFit.LineColor(colors[i%len(colors)]))
        if is_sr:
            frame_tmp.SetTitle("0-tag SIGNAL REGION [%d events]" % (datasets_all[i].sumEntries()))
            file_name = "signal_region%d.pdf" % (i-len(datasets))
        else:
            frame_tmp.SetTitle("%d-tag m_{#gamma#gamma} #in (%.2f,%.2f) [%d events]" % (args.btag, xlo, xhi, datasets_all[i].sumEntries()))
            file_name = "side_band%d.pdf" % i
        chi2 = frame_tmp.chiSquare(ndof)
        frame_tmp.Draw()
        txt = r.TPaveLabel(*(LABEL_POS+("#chi^{2}/ndof = %0.2f"%chi2,"NDC")))
        txt.SetFillStyle(0)
        txt.SetBorderSize(0)
        txt.SetTextSize(0.25)
        txt.Draw()
        up((not args.no_interact) and (not args.fastforward))
        r.gROOT.FindObject("c1").SaveAs(os.path.join(args.out,file_name))

    frame.Draw()
    up(not args.no_interact)
    r.gROOT.FindObject("c1").SaveAs(os.path.join(args.out,"sideband_pdfs.pdf"))

    fit_vals = np.array(fit_vals)
    fit_errs = np.array(fit_errs)
    fit_vals_sr = np.array(fit_vals_sr)
    fit_errs_sr = np.array(fit_errs_sr)

    # make a dictionry with all the parameter fits and errors
    # to be saved later
    param_results = {}
    for i,param_name in enumerate(param_names):
        param_results[param_name] = fit_vals[:,i]
        param_results[param_name+"_err"] = fit_errs[:,i]

    
    x_vals = arr.array('d', x_vals)
    x_errs = arr.array('d', x_errs)
    x_vals_sr = arr.array('d', x_vals_sr)
    x_errs_sr = arr.array('d', x_errs_sr)

    # for each parameter, plot the best-fit parameters vs myy slice
    r.gStyle.SetOptFit(True)
    for i,param_name in enumerate(param_names):
        p_vals = arr.array('d', fit_vals[:,i])
        p_errs = arr.array('d', fit_errs[:,i])

        g = r.TMultiGraph()

        g_sb = r.TGraphErrors(len(x_vals),x_vals,p_vals,x_errs,p_errs)

        # do a linear fit to the shape params vs. myy and
        # add the fit model and covariance to the results to be saved.
        fit_result = g_sb.Fit("pol1", "S")
        cov = fit_result.GetCovarianceMatrix()
        param_results[param_name+"_model"] = [fit_result.Parameter(0), fit_result.Parameter(1)]
        param_results[param_name+"_cov"] = [[cov[0][0], cov[0][1]], [cov[1][0], cov[1][1]]]

        g.Add(g_sb)

        if args.sr:
            # also plot any signal region points, if requested
            p_vals_sr = arr.array('d', fit_vals_sr[:,i])
            p_errs_sr = arr.array('d', fit_errs_sr[:,i])
            g_sr = r.TGraphErrors(len(x_vals_sr),x_vals_sr,p_vals_sr,x_errs_sr,p_errs_sr)
            g_sr.SetLineColor(r.kRed)
            g.Add(g_sr)

        g.SetTitle("%d-tag sideband fits;m_{#gamma#gamma} window center [GeV];best-fit %s"%(args.btag, param_name))
        g.Draw("ALP")
        g.GetYaxis().SetTitleSize(0.05)
        g.GetYaxis().SetTitleOffset(0.78)
        g.GetXaxis().SetTitleSize(0.047)
        g.GetXaxis().SetTitleOffset(0.85)
        g_sb.Draw("LP")

        up(not args.no_interact)
        file_name = "%dtag_%s_%dbins.pdf"%(args.btag,param_name,len(ds_boundaries))
        r.gROOT.FindObject("c1").SaveAs(os.path.join(args.out,file_name))

        np.savez(os.path.join(args.out, 'params.npz'),
             boundaries=ds_boundaries,
             boundary_averages=boundary_averages,
             **param_results
        )
