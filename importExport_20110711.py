#
#  importExport_20110711.py
#  
#
#  Created by Gorm Bruun Andresen on 11/07/2011.
#  Copyright (c) 2011 University of Aarhus. All rights reserved.
#
import numpy
from pylab import *
from scipy.stats.mstats import mquantiles
from Database_v1 import get_data_countries
from scipy.optimize import brentq
from scipy.optimize import fmin
from scipy.optimize import leastsq
from scipy.interpolate import Rbf
import os

#Morten code:
#from MortenStorage import get_policy_2_storage

single_column_width = 3.425
double_column_width = 2*3.425+0.236

## Works from anywhere  if you setup an ssh tunnel first: ssh -L5432:localhost:5432 gorm@pepsi.imf.au.dk
# t, L, GW, GS, datetime_offset, datalabels = get_data_countries(localhost=True);
# t, l, Gw, Gs, datetime_offset, datalabels = get_data_countries(schema='norm_agg_avg_1hour_pdata_caps_eu2020',localhost=True);

def plot_storage_vs_alpha_w(L,GW,GS,gamma=1.,alpha_w=linspace(0,1,11),path='./figures/'):

	L, GW, GS = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.
	weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

	l = weighed_sum(L)
	Gw = weighed_sum(GW)	
	Gs = weighed_sum(GS)
	
	E_H = zeros(len(alpha_w))
	for i in arange(len(alpha_w)):
		mismatch = gamma*(alpha_w[i]*Gw + (1.-alpha_w[i])*Gs) - l
		E_H[i] = get_policy_2_storage(mismatch)[1]
		
	figure(1); clf()
	
	plot(alpha_w,E_H,'k-')
	
	### Save figure:										
	save_figure('plot_storage_vs_alpha_w.pdf',gcf(),path=path)
	
	
	

def get_optimal_mix_storage(L, GW, GS, gamma=1., returnall=False):

	L, GW, GS = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.
	weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

	l = weighed_sum(L)
	Gw = weighed_sum(GW)	
	Gs = weighed_sum(GS)

	mismatch = lambda alpha_w: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l
	E_H = lambda alpha_w: get_policy_2_storage(mismatch(alpha_w))[1]
	
	alpha_w_opt = fmin(E_H,0.5,disp=True)
	
	if returnall:
	
		balancing_fixed_E_H = lambda alpha_w: sum(get_positive(-get_policy_2_storage(mismatch(alpha_w),storage_capacity = E_H(alpha_w_opt))[0])) - .01*sum(l)
	
		if balancing_fixed_E_H(0)>0:
			if balancing_fixed_E_H(1)>0:
				alpha_w_opt_1p_interval = array([brentq(balancing_fixed_E_H, 0, alpha_w_opt),brentq(balancing_fixed_E_H, alpha_w_opt, 1)])
			else:
				alpha_w_opt_1p_interval = array([brentq(balancing_fixed_E_H, 0, alpha_w_opt),1.])
		else:	
			if balancing_fixed_E_H(1)>0:
				alpha_w_opt_1p_interval = array([0.,brentq(balancing_fixed_E_H, alpha_w_opt, 1)])
			else:
				alpha_w_opt_1p_interval = array([0.,1.])
	
		print balancing_fixed_E_H(0), balancing_fixed_E_H(1)
		print alpha_w_opt, alpha_w_opt_1p_interval
	
	
		#Returns: alpha_w_opt, alpha_w_opt_1p_interval
		return alpha_w_opt, alpha_w_opt_1p_interval, E_H(alpha_w_opt)
	else:
		return E_H(alpha_w_opt), alpha_w_opt


def get_optimal_mix_balancing(L, GW, GS, gamma=1., returnall=False, normalized=True):

	L, GW, GS = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.
	weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

	l = weighed_sum(L)
	Gw = weighed_sum(GW)	
	Gs = weighed_sum(GS)
	
	mismatch = lambda alpha_w: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l
	res_load_sum = lambda alpha_w: sum(get_positive(-mismatch(alpha_w)))
	
	alpha_w_opt = fmin(res_load_sum,0.5,disp=False)
	res_load_sum_1p_interval = lambda alpha_w: res_load_sum(alpha_w)-(res_load_sum(alpha_w_opt)+.01*sum(l))
	
	alpha_w_opt_1p_interval = array([brentq(res_load_sum_1p_interval, 0, alpha_w_opt),brentq(res_load_sum_1p_interval, alpha_w_opt, 1)])
	
	if normalized:
		mismatch_opt = mismatch(alpha_w_opt)
	else:
		mismatch_opt = mismatch(alpha_w_opt)*mean(sum(L,axis=0))
	res_load_sum_opt = sum(get_positive(-mismatch_opt))
	
	if returnall:
		#Returns: alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt
		return alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt
	else:
		return alpha_w_opt


def get_mismatch_X(L, GW, GS, gamma=1.,normalized=True):

	mismatch_X = []
	for i in arange(len(L)):
		alpha_w_X, alpha_w_interval_X, res_load_sum_X, mismatch_X_ = get_optimal_mix_balancing(L[i], GW[i], GS[i], gamma, returnall=True, normalized=normalized)
		
		mismatch_X.append(mismatch_X_)
	
	return array(mismatch_X)
		
	
def get_mismatch_plus(mismatch_X):

	mismatch_plus = sum(mismatch_X*(mismatch_X>0),axis=0)

	return mismatch_plus
	
def get_mismatch_minus(mismatch_X):

	mismatch_minus = sum(mismatch_X*(mismatch_X<0),axis=0)

	return mismatch_minus


###
# plot_exchange_statistics_vs_gamma(L, GW, GS,numpy.load('country_codes.npy')['ISO'])
# 	
def plot_exchange_statistics_vs_gamma(L, GW, GS, country_labels, gamma=linspace(1,2*1.05,31), path='./figures/'):	

	results_Av = []
	for i in arange(len(gamma)):
		print i,
		sys.stdout.flush()
		
		results = plot_exchange_statistics(L, GW, GS, country_labels, gamma=gamma[i], returnall=True)
		results_Av.append(results[find(results['country_code']=='Av.')])
	
	results_Av = concatenate(results_Av)
	
	figure(1); clf()
	
	## Export
	plot(gamma,results_Av['export_last']*100/results_Av['excess'],'g--',lw=1.5)
	plot(gamma,results_Av['export_first']*100/results_Av['excess'],'g--',lw=1.5)
	pp_export = plot(gamma,results_Av['export_equal']*100/results_Av['excess'],'g-',lw=1.5)
	
	## Import
	plot(gamma,results_Av['import_last']*100/results_Av['shortage'],'r--',lw=1.5)
	plot(gamma,results_Av['import_first']*100/results_Av['shortage'],'r--',lw=1.5)
	pp_import = plot(gamma,results_Av['import_equal']*100/results_Av['shortage'],'r-',lw=1.5)
	
	axhline(100,color='k',ls='dotted')
	
	### Format figure:
	gcf().set_size_inches([single_column_width,2.8])
	#[left, bottom, width, height]
	gca().set_position([.16, .135, .835, .86], which='both')
	matplotlib.rcParams['font.size'] = 10

	axis(xmin=1,xmax=amax(gamma),ymin=0,ymax=125)
	xticks(arange(1,amax(gamma),.25))
	yticks(arange(0,105,20))
	xlabel(r'Av. RES power generation ($\gamma$) [av.l.h.]')
	ylabel(r'Av. exchange probability [%]')
	
	### Legend:
	pp = (pp_import,pp_export)
	pp_txtlabels = (r'Import',r'Export')
		
	leg = legend(pp,pp_txtlabels,loc='upper left',ncol=len(pp));
	ltext  = leg.get_texts();
	setp(ltext, fontsize='small')    # the legend text fontsize
	
	### Save figure:										
	save_figure('plot_exchange_statistics_vs_gamma.pdf',gcf(),path=path)
	
	### Reset plot default options.
	matplotlib.rcdefaults()
	
###
# plot_exchange_statistics(L, GW, GS,numpy.load('country_codes.npy')['ISO'])
# 
def plot_exchange_statistics(L, GW, GS, country_labels, gamma=1., path='./figures/',returnall=False):

	### !!!This function makes use of real units throughout most calculations!!!
	mismatch_X = get_mismatch_X(L, GW, GS, gamma, normalized=False)
	
	mismatch_plus = get_mismatch_plus(mismatch_X)
	mismatch_minus = get_mismatch_minus(mismatch_X)
	
	dtype = [('country_code','S4'),
			 ('excess','float'),
			 ('shortage','float'),
			 ('export_last','float'),
			 ('import_last','float'),
			 ('export_first','float'),
			 ('import_first','float'),
			 ('export_equal','float'),
			 ('import_equal','float')]
	results=array([zeros(len(L))],dtype=dtype,ndmin=2)
	results=results[0]
	
	results['country_code'] = country_labels
	for i in arange(len(L)):
	
		###Total need
		results['excess'][i] = mean(get_positive(mismatch_X[i]))/mean(L[i])
		results['shortage'][i] = mean(get_positive(-mismatch_X[i]))/mean(L[i])
	
		###Last in row
		i_EU_no_X = find(arange(len(L))!=i)
		mismatch_EU_no_X = sum(mismatch_X[i_EU_no_X],axis=0)
	
		results['export_last'][i] = sum(amin([abs(mismatch_X[i]),abs(mismatch_EU_no_X)],axis=0)*(mismatch_EU_no_X<0)*(mismatch_X[i]>0))/len(mismatch_X[i])/mean(L[i])
		results['import_last'][i] = sum(amin([abs(mismatch_X[i]),abs(mismatch_EU_no_X)],axis=0)*(mismatch_EU_no_X>0)*(mismatch_X[i]<0))/len(mismatch_X[i])/mean(L[i])
		
		###First in row
		results['export_first'][i] = sum(amin([abs(mismatch_X[i]),abs(mismatch_minus)],axis=0)*(mismatch_minus<0)*(mismatch_X[i]>0))/len(mismatch_X[i])/mean(L[i])
		results['import_first'][i] = sum(amin([abs(mismatch_X[i]),abs(mismatch_plus)],axis=0)*(mismatch_plus>0)*(mismatch_X[i]<0))/len(mismatch_X[i])/mean(L[i])
		
		###Equal share
		results['export_equal'][i] = sum(abs(mismatch_X[i])*amin([abs(mismatch_minus/(mismatch_plus+1e-6)),ones(L[0].shape)],axis=0)*(mismatch_X[i]>0))/len(mismatch_X[i])/mean(L[i])
		results['import_equal'][i] = sum(abs(mismatch_X[i])*amin([abs(mismatch_plus/(mismatch_minus+1e-6)),ones(L[0].shape)],axis=0)*(mismatch_X[i]<0))/len(mismatch_X[i])/mean(L[i])
	
	### Average values:
	results_Av = array([zeros(1)],dtype=dtype,ndmin=2); results_Av = results_Av[0]
	results_Av['country_code'] = 'Av.'
	results_Av['excess']       = [sum(results['excess']*sum(L,axis=1))/sum(L)]
	results_Av['shortage']     = [sum(results['shortage']*sum(L,axis=1))/sum(L)]
	results_Av['export_last']  = [sum(results['export_last']*sum(L,axis=1))/sum(L)]
	results_Av['import_last']  = [sum(results['import_last']*sum(L,axis=1))/sum(L)]
	results_Av['export_first'] = [sum(results['export_first']*sum(L,axis=1))/sum(L)]
	results_Av['import_first'] = [sum(results['import_first']*sum(L,axis=1))/sum(L)]
	results_Av['export_equal'] = [sum(results['export_equal']*sum(L,axis=1))/sum(L)]
	results_Av['import_equal'] = [sum(results['import_equal']*sum(L,axis=1))/sum(L)]
				
	### EU values:
	results_Agg = array([zeros(1)],dtype=dtype,ndmin=2); results_Agg = results_Agg[0]
	results_Agg['country_code'] = 'Agg.'
	results_Agg['excess']       = mean(get_positive(sum(mismatch_X,axis=0)))/mean(sum(L,axis=0))
	#results_Agg['excess']       = mean(mismatch_plus)/mean(sum(L,axis=0))
	#results_Agg['export_equal'] = mean(amin([abs(mismatch_plus),abs(mismatch_minus)],axis=0))/mean(sum(L,axis=0))
	
	results_Agg['shortage']     = mean(get_positive(sum(-mismatch_X,axis=0)))/mean(sum(L,axis=0))
		
	### Concatenate:
	results = concatenate([results_Agg,results_Av,results])	
	order = concatenate([[0,1],array([8,20,12,14,3,2,26,16,9,25,5,0,22,13,24,17,4,1,6,19,21,11,15,7,23,18,10])+2]) #Dominik's ordering according to lattitude.	
	x_range = concatenate([[0.,1.], 2.5+arange(len(L))])	
	
	if returnall:
		return results
		
	else:
				
		### Plot options:
		matplotlib.rcParams['xtick.direction'] = 'out'	
		matplotlib.rcParams['font.size'] = 10
		majorFormatter = FormatStrFormatter('%.2f')
				
		figure(1); clf();
		ax = gca()	#Get an axis handle
		
		### Plot data:
		bar(x_range, results['excess'][order], align='center', width=.75, color=(1,1,1), ls='dotted')	
		bar(x_range, results['export_first'][order], align='center', width=.75, color=(.8,.8,.8))	
		bar(x_range, results['export_equal'][order], align='center', width=.75, color=(.6,.6,.6))	
		bar(x_range, results['export_last'][order], align='center', width=.75, color=(.4,.4,.4))		
		
		axhline(results_Av['export_equal'],color='k',ls='--',lw=1.5)
		
		### Format figure:					
		gcf().set_size_inches([double_column_width,3])
		axis(ymin=0,ymax=0.41,xmin=amin(x_range)-.875,xmax=amax(x_range)+.875)
		xticks(x_range,results['country_code'][order],rotation=60,ha='right',va='top')
		yticks(arange(0,.41,.05))
		ylabel(r'Cumulated export [av.h.l.]')
				
		#[left, bottom, width, height]
		#pos = gca().get_position()
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		ax.spines['top'].set_color('none')
		ax.spines['right'].set_color('none')
		ax.yaxis.set_major_formatter(majorFormatter)
		ax.set_position([.085, .14, .91, .85], which='both')																								
														
		### Place legend:
		pp_total = Rectangle((0, 0), 1, 1, fill=False, ls='dotted')
		pp_first = Rectangle((0, 0), 1, 1, color=(.8,.8,.8))
		pp_equal = Rectangle((0, 0), 1, 1, color=(.6,.6,.6))
		pp_last = Rectangle((0, 0), 1, 1, color=(.4,.4,.4))
		
		pp = (pp_last,pp_equal,pp_first,pp_total)
		pp_txtlabels = (r'Last in row',r'Equal share',r'First in row',r'Excess')
		
		leg = legend(pp,pp_txtlabels,loc='upper left',ncol=len(pp));
		ltext  = leg.get_texts();
		setp(ltext, fontsize='small')    # the legend text fontsize
																																																																																																																																																																																																																																																										
		### Save figure:										
		save_figure('plot_export_statistics_gamma_{0:.0f}.pdf'.format(100*gamma),gcf(),path=path)
		
		figure(2); clf();
		ax = gca()	#Get an axis handle
		
		### Plot data:
		bar(x_range, results['shortage'][order], align='center', width=.75, color=(1,1,1), ls='dotted')	
		bar(x_range, results['import_first'][order], align='center', width=.75, color=(.8,.8,.8))	
		bar(x_range, results['import_equal'][order], align='center', width=.75, color=(.6,.6,.6))	
		bar(x_range, results['import_last'][order], align='center', width=.75, color=(.4,.4,.4))		
		
		axhline(results_Av['import_equal'],color='k',ls='--',lw=1.5)
		
		### Format figure:					
		gcf().set_size_inches([double_column_width,3])
		axis(ymin=0,ymax=.41,xmin=amin(x_range)-.875,xmax=amax(x_range)+.875)
		xticks(x_range,results['country_code'][order],rotation=60,ha='right',va='top')
		yticks(arange(0,.41,.05))
		ylabel(r'Cumulated import [av.h.l.]')
				
		#[left, bottom, width, height]
		#pos = gca().get_position()
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		ax.spines['top'].set_color('none')
		ax.spines['right'].set_color('none')
		ax.yaxis.set_major_formatter(majorFormatter)
		ax.set_position([.085, .14, .91, .85], which='both')																								
														
		### Place legend:
		pp_txtlabels = (r'Last in row',r'Equal share',r'First in row',r'Shortage')
		
		leg = legend(pp,pp_txtlabels,loc='upper left',ncol=len(pp));
		ltext  = leg.get_texts();
		setp(ltext, fontsize='small')    # the legend text fontsize
																																																																																																																																																																																																																																																										
		### Save figure:										
		save_figure('plot_import_statistics_gamma_{0:.0f}.pdf'.format(100*gamma),gcf(),path=path)	
		
		######## Normalized fractions [%]
		figure(3); clf();
		ax = gca()	#Get an axis handle
		
		### Plot data:
		bar(x_range, results['excess'][order]*100/results['excess'][order], align='center', width=.75, color=(1,1,1), ls='dotted')	
		bar(x_range, results['export_first'][order]*100/results['excess'][order], align='center', width=.75, color=(.8,.8,.8))	
		bar(x_range, results['export_equal'][order]*100/results['excess'][order], align='center', width=.75, color=(.6,.6,.6))	
		bar(x_range, results['export_last'][order]*100/results['excess'][order], align='center', width=.75, color=(.4,.4,.4))		
		
		axhline(results_Av['export_equal']*100/results_Av['excess'],color='k',ls='--',lw=1.5)
		
		### Format figure:					
		gcf().set_size_inches([double_column_width,3])
		axis(ymin=0,ymax=125,xmin=amin(x_range)-.875,xmax=amax(x_range)+.875)
		xticks(x_range,results['country_code'][order],rotation=60,ha='right',va='top')
		yticks(arange(0,105,20))
		ylabel(r'Cumulated export [%]')
				
		#[left, bottom, width, height]
		#pos = gca().get_position()
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		ax.spines['top'].set_color('none')
		ax.spines['right'].set_color('none')
		ax.yaxis.set_major_formatter(majorFormatter)
		ax.set_position([.085, .14, .91, .85], which='both')																								
														
		### Place legend:
		pp_total = Rectangle((0, 0), 1, 1, fill=False, ls='dotted')
		pp_first = Rectangle((0, 0), 1, 1, color=(.8,.8,.8))
		pp_equal = Rectangle((0, 0), 1, 1, color=(.6,.6,.6))
		pp_last = Rectangle((0, 0), 1, 1, color=(.4,.4,.4))
		
		pp = (pp_last,pp_equal,pp_first,pp_total)
		pp_txtlabels = (r'Last in row',r'Equal share',r'First in row',r'Excess')
		
		leg = legend(pp,pp_txtlabels,loc='upper left',ncol=len(pp));
		ltext  = leg.get_texts();
		setp(ltext, fontsize='small')    # the legend text fontsize
																																																																																																																																																																																																																																																										
		### Save figure:										
		save_figure('plot_export_statistics_percent_gamma_{0:.0f}.pdf'.format(100*gamma),gcf(),path=path)
		
		figure(4); clf();
		ax = gca()	#Get an axis handle
		
		### Plot data:
		bar(x_range, results['shortage'][order]*100/results['shortage'][order], align='center', width=.75, color=(1,1,1), ls='dotted')	
		bar(x_range, results['import_first'][order]*100/results['shortage'][order], align='center', width=.75, color=(.8,.8,.8))	
		bar(x_range, results['import_equal'][order]*100/results['shortage'][order], align='center', width=.75, color=(.6,.6,.6))	
		bar(x_range, results['import_last'][order]*100/results['shortage'][order], align='center', width=.75, color=(.4,.4,.4))		
		
		axhline(results_Av['import_equal']*100/results_Av['shortage'],color='k',ls='--',lw=1.5)
		
		### Format figure:					
		gcf().set_size_inches([double_column_width,3])
		axis(ymin=0,ymax=125,xmin=amin(x_range)-.875,xmax=amax(x_range)+.875)
		xticks(x_range,results['country_code'][order],rotation=60,ha='right',va='top')
		yticks(arange(0,105,20))
		ylabel(r'Cumulated import [%]')
				
		#[left, bottom, width, height]
		#pos = gca().get_position()
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		ax.spines['top'].set_color('none')
		ax.spines['right'].set_color('none')
		ax.yaxis.set_major_formatter(majorFormatter)
		ax.set_position([.085, .14, .91, .85], which='both')																								
														
		### Place legend:
		pp_txtlabels = (r'Last in row',r'Equal share',r'First in row',r'Shortage')
		
		leg = legend(pp,pp_txtlabels,loc='upper left',ncol=len(pp));
		ltext  = leg.get_texts();
		setp(ltext, fontsize='small')    # the legend text fontsize
																																																																																																																																																																																																																																																										
		### Save figure:										
		save_figure('plot_import_statistics_percent_gamma_{0:.0f}.pdf'.format(100*gamma),gcf(),path=path)
		
		#Reset plot default options.
		matplotlib.rcdefaults()
			

def plot_bar_countries(country_labels,values,x_range='EU_Agg',order='lat'):	
	
	if order=='lat':
		order = [8,20,12,14,3,2,26,16,9,25,5,0,22,13,24,17,4,1,6,19,21,11,15,7,23,18,10] #Dominik's ordering according to lattitude.
	else:	
		order = arange(len(country_labels))

	if x_range=='EU_Agg':
		x_range = concatenate([[0.,1.], 2.5+arange(len(L))])
	else:
		x_range = arange(len(country_labels))

	#Set plot options																																																																																												
	matplotlib.rcParams['xtick.direction'] = 'out'	
	matplotlib.rcParams['font.size'] = 10
	majorFormatter = FormatStrFormatter('%.2f')
	gcf().set_size_inches([double_column_width,3])
	
	ax = subplot(111)
	
	bar(x_range, values[order], align='center', width=.75, color=(.4,.4,.4), alpha=0.2)
	
	axis(ymin=0,ymax=.5,xmin=amin(x_range)-.875,xmax=amax(x_range)+.875)
	xticks(x_range,country_labels[order],rotation=60,ha='right',va='top')
	#yticks(arange(0,1.05,.25))

	ylabel(r'Wind fraction ($\alpha_w$)')
	
	#[left, bottom, width, height]
	#pos = gca().get_position()
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.spines['top'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.yaxis.set_major_formatter(majorFormatter)
	ax.set_position([.085, .14, .91, .85], which='both')

	
	
###
# plot_exchange_oportunity(L, GW, GS,numpy.load('country_codes.npy')['ISO'])
# 
# Test sample: Far distant countries
# i=[10,15,26]; plot_exchange_oportunity(L[i],GW[i], GS[i],numpy.load('country_codes.npy')['ISO'][i])
#	
# Test sample: Close countries, inland
# i=[13,22,24]; plot_exchange_oportunity(L[i],GW[i], GS[i],numpy.load('country_codes.npy')['ISO'][i])
#
#
def plot_exchange_oportunity(L, GW, GS, datalabels, gamma=1., path='./figures/'):

	mismatch_X = get_mismatch_X(L, GW, GS, gamma)
	
	mismatch_plus = get_mismatch_plus(mismatch_X)
	mismatch_minus = get_mismatch_minus(mismatch_X)

	matplotlib.rcParams['font.size'] = 10
	majorFormatter = FormatStrFormatter('%.1f')

	bin_size = .1
	for i in arange(len(L)):
		
		i_EU_no_X = find(arange(len(L))!=i)
		#alpha_w, alpha_w_interval, res_load_sum, mismatch_EU_no_X = get_optimal_mix_balancing(L[i_EU_no_X], GW[i_EU_no_X], GS[i_EU_no_X], gamma, returnall=True)
		mismatch_EU_no_X = zeros(L[0].shape)
		for k in i_EU_no_X:
		  mismatch_EU_no_X = mismatch_EU_no_X + mismatch_X[k]*sum(L[k])/sum(L[i_EU_no_X])
	
		figure(1); clf()
	
		gcf().set_dpi(300)
		gcf().set_size_inches([single_column_width,2.8])
	
		#Calculate number of bins
		x_bins = int((amax(mismatch_X[i])-amin(mismatch_X[i]))/bin_size)
		y_bins = int((amax(mismatch_EU_no_X)-amin(mismatch_EU_no_X))/bin_size)
	
		H = hexbin(mismatch_X[i],mismatch_EU_no_X,gridsize=(x_bins,y_bins),mincnt=1,cmap=matplotlib.cm.OrRd,linewidths=0)
		
		values = H.get_array()
		conf_lvl = [.9,.5,.1]
		
		levels=zeros(len(conf_lvl))
		sort_values = sort(values)
		cumsum_values = cumsum(sort_values)/sum(values)
		for j in arange(len(conf_lvl)):
			levels[j] = sort_values[argmax(cumsum_values>=(1-conf_lvl[j]))]
		
		#levels = mquantiles(values,quantiles)
		print levels
		
		paths = H.get_paths()
		xc, yc = zeros(len(values)), zeros(len(values))
		for k in arange(len(values)):
			xc[k] = mean(paths[k].vertices[0:6].transpose()[0])
			yc[k] = mean(paths[k].vertices[0:6].transpose()[1])
		
		X, Y = meshgrid(linspace(amin(xc),amax(xc),2*x_bins),linspace(amin(yc),amax(yc),2*y_bins))
		
		V = griddata(xc,yc,values,X, Y)
		V[1:-1,1:-1] = (V[1:-1,1:-1] + V[0:-2,1:-1] + V[2:,1:-1] + V[1:-1,0:-2] + V[1:-1,2:])/5.
		V[1:-1,1:-1] = (V[1:-1,1:-1] + V[0:-2,1:-1] + V[2:,1:-1] + V[1:-1,0:-2] + V[1:-1,2:])/5.
		
		#pcolor(X,Y,V)
		
		CS = contour(X,Y,V,levels,colors='k')
		
		ctxt = dict([(levels[n],'{0:.0f}%'.format(100*conf_lvl[n])) for n in arange(len(levels))])
		clabel(CS,fmt=ctxt, fontsize=8, inline=1, inline_spacing=10)
		
		
		axis(xmin=-1.5,xmax=2.3,ymin=-1.2,ymax=1.2)
		xticks(arange(-1.,2.51,1.))
		yticks(arange(-1.,1.2,.5))
		#colorbar()
		clim(0,500)
		#plot(mismatch_X[i],mismatch_EU_no_X,'k.')
	
		#Plot cross
		axhline(0,color='k',ls='-',lw=1.)
		axvline(0,color='k',ls='-',lw=1.)
	
		
	
		xlabel(r'$\Delta^{0:}$ [av.l.h.]'.format('{'+datalabels[i]+'}'))
		ylabel(r'$\Delta^{EU\backslash{} '+r'{0:}'.format(datalabels[i])+r'}$ [av.l.h.]')
		#ylabel(r'$\Delta^{EU\backslash X}$ [av.l.h.]')
	
		#[left, bottom, width, height]
		gca().set_position([.16, .135, .835, .86], which='both')
		#Text labels
		bb = gca().get_position()
		
		
		### Text labels:
		shortage_perc = sum((mismatch_X[i]<=0)*(mismatch_EU_no_X<=0))*100/len(mismatch_X[i])
		import_perc = sum((mismatch_X[i]<0)*(mismatch_EU_no_X>0))*100/len(mismatch_X[i])
		export_perc = sum((mismatch_X[i]>0)*(mismatch_EU_no_X<0))*100/len(mismatch_X[i])
		excess_perc = sum((mismatch_X[i]>=0)*(mismatch_EU_no_X>=0))*100/len(mismatch_X[i])
		
		text(bb.xmin-.14, bb.ymin-.05,'Shortage ({0:.0f}%)'.format(shortage_perc),horizontalalignment='left',verticalalignment='top',transform = gca().transAxes,weight='bold',size='small')
		text(bb.xmin-.14, bb.ymax-.1,'Import ({0:.0f}%)'.format(import_perc),horizontalalignment='left',verticalalignment='bottom',transform = gca().transAxes,weight='bold',size='small')
		text(bb.xmax-.33, bb.ymin-.05,'Export ({0:.0f}%)'.format(export_perc),horizontalalignment='left',verticalalignment='top',transform = gca().transAxes,weight='bold',size='small')
		text(bb.xmax-.33, bb.ymax-.1,'Excess ({0:.0f}%)'.format(excess_perc),horizontalalignment='left',verticalalignment='bottom',transform = gca().transAxes,weight='bold',size='small')
		
		gca().xaxis.set_major_formatter(majorFormatter)
		gca().yaxis.set_major_formatter(majorFormatter)
	
		save_figure('plot_exchange_oportunity_last_'+datalabels[i]+'_gamma_{0:.0f}.pdf'.format(gamma*100),gcf(),path=path)

		#figure(2); clf()
		
		##First in row
		#mismatch_plus_minus = (mismatch_X[i]<=0)*(mismatch_plus*(mismatch_plus>0) + mismatch_minus*(mismatch_plus<=0)) + \
		#	(mismatch_X[i]>0)*(mismatch_minus*(mismatch_minus<0) + mismatch_plus*(mismatch_minus>=0))
		#mismatch_plus_minus = mismatch_plus_minus/len(L)


		##Calculate number of bins
		#x_bins = int((amax(mismatch_X[i])-amin(mismatch_X[i]))/bin_size)
		#y_bins = int((amax(mismatch_plus_minus)-amin(mismatch_plus_minus))/bin_size)

		##Plot 2D histogram
		#hexbin(mismatch_X[i],mismatch_plus_minus,gridsize=(x_bins,y_bins),mincnt=1)
		#axis(xmin=-1.5,xmax=2.5,ymin=-1.5,ymax=1.5)
		#colorbar()
		##clim(0,500)

		##Plot cross
		#axhline(0,color='k',ls='-',lw=1.5)
		#axvline(0,color='k',ls='-',lw=1.5)

		#save_figure('plot_exchange_oportunity_1st_'+datalabels[i]+'.pdf',gcf(),path=path)

	return H

def plot_storage_vs_gamma(L, GW, GS, gamma=linspace(1.,2.*1.05,21),path='./figures/'):

	E_H_X = zeros((len(gamma),len(L)))
	E_H_Av, E_H_EU = zeros(len(gamma)), zeros(len(gamma))
	for i in arange(len(gamma)):
		for k in arange(len(L)):
			E_H_X[i,k], alpha_w  = get_optimal_mix_storage(L[k], GW[k], GS[k], gamma[i])

		E_H_Av[i] = sum(E_H_X[i]*sum(L,axis=1))/sum(L)
		E_H_EU[i], alpha_w_EU  =  get_optimal_mix_storage(L, GW, GS, gamma[i])

	matplotlib.rcParams['font.size'] = 10

	figure(1); clf()

	gcf().set_dpi(300)
	gcf().set_size_inches([single_column_width,2.8])

	pp_Av = plot(gamma,E_H_Av/365./24.,'k-')
	pp_EU = plot(gamma,E_H_EU/365./24.,'k--')
	
	xticks(arange(amin(gamma),amax(gamma)*1.05,.25))
	yticks(arange(0,.2*1.05,.05))
	axis(xmin=1.,xmax=2.*1.05,ymin=0,ymax=.2*1.05)
	
	xlabel(r'Av. RES power generation ($\gamma$) [av.l.h.]')
	ylabel(r'Storage energy cap. ($E_H$) [av.l.yr.]')
	
	place_legend((pp_Av,pp_EU),(r'Av.',r'Agg.'),loc='upper right')
	
	#[left, bottom, width, height]
	gca().set_position([.17, .135, .825, .86], which='both')
	
	save_figure('plot_storage_vs_gamma.pdf',gcf(),path=path)
	
	
def plot_balancing_vs_gamma(L, GW, GS, gamma=linspace(1.,2.*1.05,11),path='./figures/'):

	res_load_sum_X, res_load_q90_X, res_load_q99_X = zeros((len(gamma),len(L))), zeros((len(gamma),len(L))), zeros((len(gamma),len(L)))
	L_quantile =  zeros(len(L))
	res_load_sum_Av, res_load_q90_Av, res_load_q99_Av = zeros(len(gamma)), zeros(len(gamma)), zeros(len(gamma))
	res_load_sum_EU, res_load_q90_EU, res_load_q99_EU = zeros(len(gamma)), zeros(len(gamma)), zeros(len(gamma))
	for i in arange(len(gamma)):
		for k in arange(len(L)):
			alpha_w, alpha_w_1p_interval, res_load_sum_X[i,k], mismatch = get_optimal_mix_balancing(L[k], GW[k], GS[k], gamma[i], returnall=True)
			res_load_q90_X[i,k], res_load_q99_X[i,k] = float(mquantiles(-mismatch,.90)), float(mquantiles(-mismatch,.99))
			
			if i==0:
				L_quantile[k] = float(mquantiles(L[k],.99)[0]/mean(L[k]))
	
		#Average
		res_load_sum_Av[i] = sum(res_load_sum_X[i]*sum(L,axis=1))/sum(L)
		res_load_q90_Av[i] = sum(res_load_q90_X[i]*sum(L,axis=1))/sum(L)
		res_load_q99_Av[i] = sum(res_load_q99_X[i]*sum(L,axis=1))/sum(L)
	
		#Aggregated
		alpha_w_EU, alpha_w_1p_interval_EU, res_load_sum_EU[i], mismatch =  get_optimal_mix_balancing(L, GW, GS, gamma[i], returnall=True)		
		res_load_q90_EU[i], res_load_q99_EU[i] = float(mquantiles(-mismatch,.90)), float(mquantiles(-mismatch,.99))
	
	L_quantile_Av = [sum(L_quantile*sum(L,axis=1))/sum(L)]																								
	L_quantile_EU = mquantiles(sum(L,axis=0),.99)[0]/mean(sum(L,axis=0))

	matplotlib.rcParams['font.size'] = 10

	figure(1); clf()

	gcf().set_dpi(300)
	gcf().set_size_inches([single_column_width,2.8])

	pp_EU = plot(gamma,res_load_sum_EU,'k-')
	pp_Av = plot(gamma,res_load_sum_Av,'k--')
	
	xticks(arange(1.,2.*1.05,.25))
	yticks(arange(0,.25*1.05,.05))
	axis(xmin=1.,xmax=2.*1.05,ymin=0,ymax=.25*1.05)
	
	xlabel(r'Av. RES power generation ($\gamma$) [av.l.h.]')
	ylabel(r'Av. residual load ($\langleR\rangle_t$) [av.l.h.]')
	
	place_legend((pp_Av,pp_EU),(r'Av.',r'Agg.'),loc='upper right')
	
	#[left, bottom, width, height]
	gca().set_position([.17, .135, .825, .86], which='both')
	
	save_figure('plot_balancing_vs_gamma.pdf',gcf(),path=path)
	
	figure(2); clf() 
	matplotlib.rcdefaults()
	gcf().set_dpi(300)
	gcf().set_size_inches([single_column_width,2.8])

	pp_EU = plot(gamma,res_load_q99_EU,'k-')  #Shared RES and balancing capacity.
	pp_Av = plot(gamma,res_load_q99_Av,'k--')   #No sharing.
	pp_L = axhline(L_quantile_Av,color='k',ls=':',lw=1.5)
	
	xticks(arange(1.,2.*1.05,.25))
	yticks(arange(0,1.5*1.05,.25))
	axis(xmin=1.,xmax=2.*1.05,ymin=0,ymax=1.5*1.05)
	
	xlabel(r'Av. RES power generation ($\gamma$) [av.l.h.]')
	ylabel(r'Balancing power cap. [av.l.h.]')
	
	leg = legend((pp_L,pp_Av,pp_EU),(r'$Q_{99\%}(L)$',r'$Q_{99\%}(R_{Av.})$',r'$Q_{99\%}(R_{Agg.})$'),loc='lower left',ncol=3,columnspacing=.1,handlelength=2,handletextpad=0);
	ltext  = leg.get_texts();
	setp(ltext, fontsize='small')    # the legend text fontsize
	
	#[left, bottom, width, height]
	gca().set_position([.17, .135, .825, .86], which='both')
	
	save_figure('plot_quantiles_vs_gamma.pdf',gcf(),path=path)

###
# plot_optimal_mix_storage_bar(L, GW, GS,numpy.load('country_codes.npy')['ISO'])
# plot_optimal_mix_storage_bar(L, GW, GS,numpy.load('country_codes.npy')['ISO'], gamma=1.5)
#
def plot_optimal_mix_storage_bar(L, GW, GS, datalabels, gamma=1.,path='./figures/'):

	order = [8,20,12,14,3,2,26,16,9,25,5,0,22,13,24,17,4,1,6,19,21,11,15,7,23,18,10] #Latitude
	order = get_order_RES_quality(GW,GS,datalabels) #wind cap. factor/solar cap. factor

	alpha_w, alpha_w_1p_interval, E_H = zeros(len(L)), list(zeros(len(L))), zeros(len(L))
	for i in arange(len(L)):
		alpha_w[i], alpha_w_1p_interval[i], E_H[i] = get_optimal_mix_storage(L[i], GW[i], GS[i], gamma, returnall=True)
	
	return alpha_w, alpha_w_1p_interval, order
	alpha_w_1p_interval = array(alpha_w_1p_interval)
		
	x_range = concatenate([[0.,1.], 2.5+arange(len(L))]);	
	datalabels = concatenate([['Agg.','Av.'],datalabels[order]])
	
	alpha_w_EU, alpha_w_1p_interval_EU, E_H_EU =  get_optimal_mix_storage(L, GW, GS, gamma, returnall=True)		
	
	alpha_w_Av = [sum(alpha_w*sum(L,axis=1))/sum(L)]						
	alpha_w_1p_interval_Av = [sum(alpha_w_1p_interval.transpose()[0]*sum(L,axis=1))/sum(L),sum(alpha_w_1p_interval.transpose()[1]*sum(L,axis=1))/sum(L)]
	E_H_Av = [sum(E_H*sum(L,axis=1))/sum(L)]

	#Set plot options																																																																																												
	matplotlib.rcParams['xtick.direction'] = 'out'	
	matplotlib.rcParams['font.size'] = 10
	majorFormatter = FormatStrFormatter('%.2f')
	
	#Make figure																																					
	figure(1); clf()
	
	gcf().set_dpi(300)
	gcf().set_size_inches([double_column_width,3])
	
	ax = subplot(111)
	
	alpha_w_1p_interval = concatenate([[alpha_w_1p_interval_EU], [alpha_w_1p_interval_Av], alpha_w_1p_interval[order]]).transpose()
	alpha_w = concatenate([alpha_w_EU,alpha_w_Av,alpha_w[order]])
	alpha_w_1p_interval_ = [alpha_w - alpha_w_1p_interval[0], alpha_w_1p_interval[1] - alpha_w]

	bar(x_range, alpha_w, yerr=alpha_w_1p_interval_, ecolor='k', error_kw={'elinewidth':1.}, align='center', width=.75, color=(.4,.4,.4))
	axhline(alpha_w_Av,color='k',ls='--',lw=1.5)
	
	axis(ymin=0,ymax=1.05,xmin=amin(x_range)-.875,xmax=amax(x_range)+.875)
	xticks(x_range,datalabels,rotation=60,ha='right',va='top')
	yticks(arange(0,1.05,.25))

	ylabel(r'Wind fraction ($\alpha_w$)')
	
	#[left, bottom, width, height]
	#pos = gca().get_position()
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.spines['top'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.yaxis.set_major_formatter(majorFormatter)
	ax.set_position([.085, .14, .91, .85], which='both')
	
	save_figure('plot_optimal_mix_storage_bar_gamma{0:.0f}_cap_factor.pdf'.format(gamma*100),gcf(),path=path)

	#Make figure																																					
	figure(2); clf()
	
	gcf().set_dpi(300)
	gcf().set_size_inches([double_column_width,3])
	
	ax = subplot(111)
	
	E_H = concatenate([[E_H_EU],E_H_Av,E_H[order]])

	bar(x_range, E_H/365./24., align='center', width=.75, color=(.4,.4,.4))
	axhline(array(E_H_Av)/365./24.,color='k',ls='--',lw=1.5)
	
	axis(ymin=0,ymax=.3*1.05,xmin=amin(x_range)-.875,xmax=amax(x_range)+.875)
	xticks(x_range,datalabels,rotation=60,ha='right',va='top')
	yticks(arange(0,.3*1.05,.05))

	ylabel(r'Storage energy capacity ($E_H$) [av.l.yr.]')
	
	#[left, bottom, width, height]
	#pos = gca().get_position()
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.spines['top'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.yaxis.set_major_formatter(majorFormatter)
	ax.set_position([.085, .14, .91, .85], which='both')
	
	save_figure('plot_minimal_storage_bar_gamma{0:.0f}_cap_factor.pdf'.format(gamma*100),gcf(),path=path)

	#Reset plot default options.
	matplotlib.rcdefaults()



###
# plot_optimal_mix_bar(L, GW, GS,numpy.load('country_codes.npy')['ISO'])
#
#
def plot_optimal_mix_bar(L, GW, GS, datalabels, gamma=1.,path='./figures/'):

	order = [8,20,12,14,3,2,26,16,9,25,5,0,22,13,24,17,4,1,6,19,21,11,15,7,23,18,10]													

	alpha_w, alpha_w_1p_interval, res_load_sum, res_load_q90, res_load_q99, L_quantile = zeros(len(L)), list(zeros(len(L))), zeros(len(L)), zeros(len(L)), zeros(len(L)), zeros(len(L))
	for i in arange(len(L)):
		alpha_w[i], alpha_w_1p_interval[i], res_load_sum[i], mismatch = get_optimal_mix_balancing(L[i], GW[i], GS[i], gamma, returnall=True)
		res_load_q90[i], res_load_q99[i] = mquantiles(-mismatch,.90), mquantiles(-mismatch,.99)
		L_quantile[i] = mquantiles(L[i],.99)[0]/mean(L[i])
	alpha_w_1p_interval = array(alpha_w_1p_interval)
		
	x_range = concatenate([[0.,1.], 2.5+arange(len(L))]);	
	datalabels = concatenate([['Agg.','Av.'],datalabels[order]])
	
	alpha_w_EU, alpha_w_1p_interval_EU, res_load_sum_EU, mismatch =  get_optimal_mix_balancing(L, GW, GS, gamma, returnall=True)		
	res_load_q90_EU, res_load_q99_EU = mquantiles(-mismatch,.90), mquantiles(-mismatch,.99)
	L_quantile_EU = mquantiles(sum(L,axis=0),.99)[0]/mean(sum(L,axis=0))

	alpha_w_Av = [sum(alpha_w*sum(L,axis=1))/sum(L)]						
	alpha_w_1p_interval_Av = [sum(alpha_w_1p_interval.transpose()[0]*sum(L,axis=1))/sum(L),sum(alpha_w_1p_interval.transpose()[1]*sum(L,axis=1))/sum(L)]
	res_load_sum_Av = [sum(res_load_sum*sum(L,axis=1))/sum(L)]
	L_quantile_Av = [sum(L_quantile*sum(L,axis=1))/sum(L)]																								
	res_load_q90_Av = [sum(res_load_q90*sum(L,axis=1))/sum(L)]
	res_load_q99_Av = [sum(res_load_q99*sum(L,axis=1))/sum(L)]
	
	#Set plot options																																																																																												
	matplotlib.rcParams['xtick.direction'] = 'out'	
	matplotlib.rcParams['font.size'] = 10
	majorFormatter = FormatStrFormatter('%.2f')
	
																																						
	figure(1); clf()
	
	gcf().set_dpi(300)
	gcf().set_size_inches([double_column_width,3])
	
	ax = subplot(111)
	
	alpha_w_1p_interval = concatenate([[alpha_w_1p_interval_EU], [alpha_w_1p_interval_Av], alpha_w_1p_interval[order]]).transpose()
	alpha_w = concatenate([alpha_w_EU,alpha_w_Av,alpha_w[order]])
	alpha_w_1p_interval_ = [alpha_w - alpha_w_1p_interval[0], alpha_w_1p_interval[1] - alpha_w]

	bar(x_range, alpha_w, yerr=alpha_w_1p_interval_, ecolor='k', error_kw={'elinewidth':1.}, align='center', width=.75, color=(.4,.4,.4))
	axhline(alpha_w_Av,color='k',ls='--',lw=1.5)
	
	axis(ymin=0,ymax=1.05,xmin=amin(x_range)-.875,xmax=amax(x_range)+.875)
	xticks(x_range,datalabels,rotation=60,ha='right',va='top')
	yticks(arange(0,1.05,.25))

	ylabel(r'Wind fraction ($\alpha_w$)')
	
	#[left, bottom, width, height]
	#pos = gca().get_position()
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.spines['top'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.yaxis.set_major_formatter(majorFormatter)
	ax.set_position([.085, .14, .91, .85], which='both')
	
	save_figure('plot_optimal_mix_res_load_bar_gamma_{0:.2f}.pdf'.format(100*gamma),gcf(),path=path)
	
	figure(2); clf()
	
	gcf().set_dpi(300)
	gcf().set_size_inches([double_column_width,3])
	
	ax = subplot(111)
	
	bar(x_range,concatenate([[L_quantile_EU],L_quantile_Av,L_quantile[order]]), align='center', width=.75, fill=False, ls='dotted')
	bar(x_range,concatenate([res_load_q99_EU,res_load_q99_Av,res_load_q99[order]]), align='center', width=.75, color=(.8,.8,.8))
	bar(x_range,concatenate([res_load_q90_EU,res_load_q90_Av,res_load_q90[order]]), align='center', width=.75, color=(.6,.6,.6))
	bar(x_range,concatenate([[res_load_sum_EU],res_load_sum_Av,res_load_sum[order]])/len(L[0]), align='center', width=.75, color=(.4,.4,.4))
	axhline(array(res_load_sum_Av)/len(L[0]),color='k',ls='--',lw=1.5)

	pp_load_q99 = Rectangle((0, 0), 1, 1, fill=False, ls='dotted')
	pp_res_load_q99 = Rectangle((0, 0), 1, 1, color=(.8,.8,.8),lw=1)
	pp_res_load_q90 = Rectangle((0, 0), 1, 1, color=(.6,.6,.6))
	pp_res_load_Av = Rectangle((0, 0), 1, 1, color=(.4,.4,.4))

	#plot(x_range,concatenate([[L_quantile_EU],L_quantile_Av,L_quantile[order]]))
	#plot(x_range,concatenate([res_load_q90_EU,res_load_q90_Av,res_load_q90[order]]))
	#plot(x_range,concatenate([res_load_q99_EU,res_load_q99_Av,res_load_q99[order]]))
	
	axis(ymin=0,ymax=2.1,xmin=amin(x_range)-.875,xmax=amax(x_range)+.875)
	xticks(x_range,datalabels,rotation=60,ha='right',va='top')
	yticks(arange(0,2.1,.5))
	ylabel(r'Residual load [av.l.h.]')
	
	pp = (pp_res_load_Av,pp_res_load_q90,pp_res_load_q99,pp_load_q99)
	pp_txtlabels = (r'$\langle R \rangle_t$',r'$Q_{90\%}(R)$',r'$Q_{99\%}(R)$',r'$Q_{99\%}(L)$')
	
	leg = legend(pp,pp_txtlabels,loc='upper left',ncol=len(pp));
	ltext  = leg.get_texts();
	setp(ltext, fontsize='small')    # the legend text fontsize
	
	#place_legend((pp_res_load_Av,pp_res_load_q90,pp_res_load_q99,pp_load_q99),(r'$R_{mean}$',r'$Q_{90\%}(R)$',r'$Q_{99\%}(R)$',r'$Q_{99\%}(L)$'),loc='upper right')
	
	#[left, bottom, width, height]
	#pos = gca().get_position()
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.spines['top'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.yaxis.set_major_formatter(majorFormatter)
	ax.set_position([.085, .14, .91, .85], which='both')
	
	save_figure('plot_minimal_res_load_bar_gamma_{0:.2f}.pdf'.format(100*gamma),gcf(),path=path)

	#Reset plot default options.
	matplotlib.rcdefaults()

######
#  plot_region_time_series(t, L, GW, GS, txtlabel='EU')
#  for i in arange(27): plot_region_time_series(t, L[i], GW[i], GS[i], txtlabel=datalabels[i])
#
def plot_region_time_series(t, L, GW, GS, tau=30*24, path='./figures/', txtlabel=''):

	L, GW, GS = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.
	weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))
		
	color_wind = (0.5,0.7,1.)
	color_solar = (1.,.8,0.)

	t = t + datestr2num('1-1-2000')

	load = weighed_sum(L)
	load_slow, load_ = get_high_low_fft(load, tau)

	Gw = weighed_sum(GW)
	Gw_slow, Gw_ = get_high_low_fft(Gw, tau)

	Gs = weighed_sum(GS)
	Gs_slow, Gs_ = get_high_low_fft(Gs, tau)

	figure(1); clf()

#	plot(t,Gs,'.',color=color_solar,alpha=0.5,ms=10,markeredgewidth=0)
#	plot(t,Gw,'.',color=color_wind,alpha=0.5,ms=8,markeredgewidth=0)
	
	plot(t,Gs,'-',color=color_solar,alpha=0.5,aa=True)
	plot(t,Gw,'-',color=color_wind,alpha=.9,aa=True)
	plot(t,load,'-',color='k',alpha=0.6,aa=True)
	
	plot(t,load_slow,'w',lw=5.5)	
	pp_load = plot(t,load_slow,'k',lw=4)	

	plot(t,Gs_slow,'w',lw=6)
	pp_solar = plot(t,Gs_slow,'-',color=color_solar,lw=4)
	
	plot(t,Gw_slow,'w',lw=5.5)
	pp_wind = plot(t,Gw_slow,'-',color=color_wind,lw=4)

	t_min = datestr2num('31-12-2000')
	t_max = datestr2num('31-12-2001')
	t_range = t[find((t>=t_min) * (t<=t_max))]
	
	i_months = find(diff([num2date(t).month for t in t_range])) + 1
	
	axis(ymin=0,ymax=6.5,xmin=t_min,xmax=t_max)
	xticks(t_range[i_months],[num2date(t).strftime('%b') for t in t_range[i_months]],rotation=-45,ha='left')
	yticks(arange(7))
	ylabel('Normalized power [av.h.l.]')

	place_legend((pp_load[0],pp_wind[0],pp_solar[0]),('Load','Wind','Solar PV'),loc='upper right')

	save_figure('plot_region_load_'+txtlabel+'.pdf',path=path)

def save_country_codes():

	names = ['Austria','Belgium','Bulgaria','Bosnia and Herzegovina','Czech Republic','Switzerland','Germany','Denmark','Spain','France','Finland','Great Britain',\
		'Greece','Hungary','Italy','Ireland','Croatia','Luxembourg','Norway','Netherlands','Portugal','Poland','Romania','Sweden','Slovakia','Slovenia','Serbia']
	ISO_codes  = ['AT','BE','BG','BA','CZ','CH','DE','DK','ES','FR','FI','GB','GR','HU','IT','IE','HR','LU','NO','NL','PT','PL','RO','SE','SK','SI','RS']
	ISET_codes  = ['A','B','BG','BH','CZ','Ch','D','DK','ES','F','FIN','GB','GR','H','I','IRL','Kro','Lux','N','NL','P','PL','Ro','S','SK','SLO','SRB']
	
	table = array([empty(len(names))],dtype={'names': ('name', 'ISO', 'ISET'),'formats': ('S30', 'S8', 'S8')})
	table['name'] = names
	table['ISO'] = ISO_codes
	table['ISET'] = ISET_codes
 
	numpy.save('country_codes',table[0])


###
# Utilities: #
###

def save_figure(figname='TestFigure.png', fig=gcf(), path='./', dpi=300):
	
	figure(fig.number)
	savefig(path + figname, dpi=300)
	print 'Saved figure:',path + figname
	sys.stdout.flush()
	
def place_legend(pp,txtlabels,loc='upper right'):

	leg = legend(pp,txtlabels,loc=loc);
	ltext  = leg.get_texts();
	setp(ltext, fontsize='small')    # the legend text fontsize
	
def get_high_low_fft(load, tau=24):
	"Error if len(load) not even"
	
	iseven = 2*(len(load)/2) == len(load)
	if not iseven:
		print "Warning! Uneven length array. Patching the end. (dfkl32k)"
		load = concatenate([load,[load[-1]]])
	
	time = arange(len(load));
	freq = fftfreq(len(time))[:len(time)/2+1]
	load_fft = rfft(load);

	sigma = tau*2;
	kernel = exp(-2*(pi*freq/(2*pi)*sigma)**2);

	load_fft = rfft(load);
	load_high_gauss = irfft(load_fft*(1-kernel)).real;
	load_low_gauss = load - load_high_gauss;
	
	if iseven:
		return load_low_gauss, load_high_gauss
	else:
		return load_low_gauss[:-1], load_high_gauss[:-1]
		
def get_positive(x):
	
	return x*(x>0)
	
def load_country_capacities(filename='Iset_country_capacities.txt'):

	return  np.loadtxt('Iset_country_capacities.txt',skiprows=1,dtype=[('country','S4'),('WP_cap','float'),('PV_cap','float')])

def get_country_capacity_factors(GW,GS,datalabels,returnall=False):

	country_capacities = load_country_capacities()
	
	if not all([datalabels[i]==country_capacities['country'][i] for i in arange(len(datalabels))]):
		print "WARNING! Country label mismatch!!! (m23ncx)"
	
	wind_factors = mean(GW, axis=1)/(country_capacities['WP_cap']*1e-3)
	wind_std = std(GW, axis=1)/(country_capacities['WP_cap']*1e-3)
	
	solar_factors = mean(GS, axis=1)/(country_capacities['PV_cap']*1e-3)
	solar_std = std(GS, axis=1)/(country_capacities['PV_cap']*1e-3)
	
	if returnall:
		return wind_factors, solar_factors, wind_std, solar_std
	else:
		return wind_factors, solar_factors
	
def get_order_RES_quality(GW,GS,datalabels):

	wind_factors, solar_factors = get_country_capacity_factors(GW,GS,datalabels)
	
	ratio = wind_factors/solar_factors
	
	return argsort(ratio)
	
def plot_aoptimal_mix_vs_RES_rank(GW,GS,datalabels):

	wind_factors, solar_factors = get_country_capacity_factors(GW,GS,datalabels)

	rank = wind_factors - solar_factors
	
	order_rank = argsort(rank)
