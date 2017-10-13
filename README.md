bbyy-sidebands
===

This is a tool to estimate the shape parameters of the bbyy 4-body mass distribution.
The data is sliced into bins of m_yy, and the specified model (novosibirsk, landau, etc), is fit to the m_yybb distribution in each slice.
The best-fit shape parameters of each slice are then plotted vs. m_yy and linear model is fit that may be interpolated to the signal region.

Usage
===

The script needs MxAOD data to run on, that has been pre-selected up to (but not including) the tight photon mass cut.
If you're on lxplus this is available in Chase's public area.

Options that you might care about:

 * `--binning`: The binning scheme to use, which defines the boundaries of the individual m_yy slices. `sr-window-1` is defined such that exactly one slice fits in the signal region, `sr-window-2` has two slices in the SR, etc. You can define your own towards the beginning of the source. See `-h` for available options.
 * `--btag`: The btag category to select and fit to (0,1,2).
 * `--data-path`: path containing one or more `.root` files w/ the skimmed MxAOD data.
 * `--out`: path to save all the plots as well as a numpy file with the best-fit parameter information.
 * `--sr`: Also fit slices in the signal region (only works when `--btag 0` is specified).

Options you probably don't care about:
 * `--model`: The BG pdf model to use (default is Novosibirsk). See `-h` for available options.
 * `--fastforward`: don't pause interactively after each fit, just loop until all the m_yy regions are done.
 * `--no-interact`: skip all the interactive prompts.
 * `--highmass`: use the high mass selection (by default does the low mass analyses)
 * `--tail`: fix the Novosibirsk tail parameter (which sometimes in finnicky in fits) to a specific value.
