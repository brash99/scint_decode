LIST OF FUNCTIONS in init.C

You may notice that most of the plot functions names may or may have
not a "1" at the end of the name.  The versions *with* the "1"
(e.g. plot_eta1) give a single plot, as opposed to the versions
*without* the "1" (e.g. plot_eta), which give a multiplot of 16
channels.



//
// RUN THIS FIRST!
//



init(runno)
// reads root file corresponding to run runno, extracts the tree with raw data,
// define shortcuts for slot 0 AT THE MOMENT, ONLY SLOT 0 CHANS 0-47 HAVE 
// SHORTCUTS!!!
// Note: in the general function section I added functions to convert slots and
// chans ADC from/to TDC. They are used in the plot_* routines, but I am not sure
// how to meaningfully use them for the following aliases



//
// PLOT_ functions
//



plot_eta1( adc_slot, adc_chan, tdc_cut, tdc_slot = -1, tdc_chan = -1, Ath_suggested=75, Aw_suggested=7)
// plot eta (calculated from histos) and its fit for a given slot and
// channel.
// Starting values for threshold and width may be suggested.
// The fit results are stored in global arrays Athresholds and Awidths

plot_eta( adc_slot, adc_chan_start, tdc_cut, tdc_slot = -1, tdc_chan_start = -1)
// plot 16 eta (calculated from histos) and its fit for a given slot
// and from a channel. The fit results are stored in global arrays Athresholds 
// and Awidths

plot_tdc( tdc_slot, tdc_chan_start, adc_cut, adc_slot = -1, adc_chan_start = -1)
// plot canvas divided in 4x4 TDC plots of 16 channels starting from
// adc_chan_start, of given slot. TDC is plotted with and without the
// following cut: adc value > adc_cut
// By default, tdc starting chan = adc starting chan; if passed a
// different value than -1 to adc_chan_start, that value will be used
// instead. May be useful during testing of ADC/TDC channels
// Analogously for adc_slot

plot_adc( adc_slot, adc_chan_start, tdc_cut, tdc_slot = -1, tdc_chan_start = -1
// plot canvas divided in 4x4 ADC plots of 16 channels starting from
// adc_chan_start, of given slot. ADC is plotted with and without the
// following cut: tdc value > tdc_cut
// By default, tdc starting chan and slot are calculated
// automatically; if a different values than -1 to
// tdc_chan_start/slot, those values will be used instead. May be
// useful during testing of ADC/TDC channels

plot_tdc_vs_adc(adc_slot, adc_chan_start, tdc_slot = -1, tdc_chan_start = -1){
// For interactive use: plot canvas divided in 4x4 TDC vs ADC plots of
// 16 channels starting from adc_chan_start, of given slot. No cuts.
// By default, tdc starting chan is calculated from adc starting chan;
// if passed a different value than -1 to tdc_chan_start, that value
// will be used instead. May be useful during testing of ADC/TDC
// channels. Fills global vector for noise freq

plot_peak1( slot, chan, value_lt=10, value_gt=100)
// plots the peak for vertical tracks for given slot/channel imposing
// the adjacent channels have counts < value_lt, and the chosen
// channel at least value_gt counts. Fills the global vectors

plot_peak( slot, chan_start, value_lt=10, value_gt=100)
// plots the peak for vertical tracks for 16 chans starting from given
// slot/channel imposing the adjacent channels have counts < value_lt,
// and the chosen channel at least value_gt counts. Fills the global vectors

plot_walk1(int adc_slot, int adc_chan, int tdc_slot=-1, int tdc_chan=-1, int adc_min=0, int adc_max=600, int deltat_min=0, int deltat_max=120, int stage=6)
// plots the walk correction of given adc(tdc) slot, channel, for adc
// in the range [adc_min, adc_max], and delta t (= leading edge -
// trailing edge) in range [deltat_min, deltat_max].
//
// The stage flag may be used to see various steps of the calculation:
// stage == 1 user wants to check raw data;
// stage == 2 user wants to see the profile histogram of stage 1;
// stage == 3 user wants to see the distribution of ADC - mean;
// stage == 4 user wants to see ADC-mean vs ADC;
// stage == 5 user wants to see (ADC-mean)/ADC vs ADC;
// stage == 6 same as stage == 5, but profile



//
// AUXiliary functions. Normally there is no reason to call them
// interactively
//



noise_freq(tdc_slot, tdc_chan, tdc_cut)
// calculate the noise frequency (in Hz) for given tdc slot, chan, cut

phe(adc_slot, adc_chan)
// Calculate the number of phe for a given adc slot, chain

gain(adc_slot, adc_chan)
// calculate gain for a given adc slot, chan

calc_adc_chan( tdc_slot, tdc_chan)
// find the adc slot corresponding to a given tdc slot/channel

calc_adc_slot( tdc_slot, tdc_chan)
// find the tdc slot corresponding to a given tdc slot/channel

calc_tdc_chan( adc_slot, adc_chan)
// find the tdc chan corresponding to a given adc slot/channel

calc_tdc_slot( adc_slot, adc_chan)
// find the tdc slot corresponding to a given adc slot/channel

fit_eta_in_slot( adc_slot=0, adc_chan_min=0, adc_chan_max=63, tdc_cut=1200, tdc_slot = -1, tdc_chan = -1)
// auxiliary function. Calculate all eta fits for slot, between
// chan_min and chan_max (inclusive), in order to populate the
// arrays. Used in save_eta()

peak_rms( slot, chan, value_lt=10, value_gt=100)
// find the RMS of the peak for a given slot/channel imposing the
// adjacent channels have counts < value_lt, and the chosen channel
// at least value_gt counts

peak_mean( slot, chan, value_lt=10, value_gt=100)
// find the position of the peak for a given slot/channel imposing the
// adjacent channels have counts < value_lt, and the chosen channel
// at least value_gt counts

save_data(abs_NINOthreshold_V, abs_PMT_HV_V, 
	       adc_slot=-1, chan_min=0, chan_max=63, 
	       tdc_slot=-1, tdc_chan=-1)
// For interactive use:
// save results (retrieved from global arrays, see "format" variable for
// details) on a file named after run number, NINO threshold and PMT
// voltage (overwrite mode).
// Default values meaning:
// save every ADC slot (slot = -1), every ADC channel (0,63). 
// TDC slot/chan are calculated from ADC ones.
