# coding=utf-8
import numpy as np
import pandas as pd
import numpy_financial as npf
import matplotlib.pyplot as plt

from DK_Collections import basic_plot_polishing, plot_exit

def est_systemcosts(kW_DC):
	"""Returns an estimate of the system costs in 2020 USD (i.e. investment costs) given the size of the PV system in
	kW DC.

	Costs:
		PV panels
		Inverter
		BOS (excluding inverter)
		Installation labor costs
		Other soft costs (tax, developer and EPC overheads, PII)

	References:
		https://www.nrel.gov/docs/fy19osti/72399.pdf   2018 benchmark prices in the US (*1.02 to adjust to 2020 dollars)
	"""
	# 2020 USD/W DC
	_default = pd.Series({
		'PV panels': 0.47,
		'inveter': 0.08,
		'BOS (w/o inverter)': 0.24,
		'installation labor': 0.16,
		'other soft costs': 0.87,
	}, name='USD 2020')*1.02

	return _default * kW_DC * 10**3


def calc_cashflows(InvCosts, Sys_kW, solar_yield, energy_price, L, ResaleValue, deg_early=0.0,
                   deg_cons=0.8,):
	"""Calculates the cashflows of a solar PV system.

	PARAMETERS:
		InvCosts            Total investment costs, paid at the beginning of the project. This includes the total
							system costs (suggest to use est_systemcosts()), and any tax credits.

		Sys_kW              The system capacity in kW DC
		solar_yield         The annual solar yield at the site.
		deg_early           Argument of  calc_annualprod()
		deg_cons            Argument of  calc_annualprod()

		energy_price        The price of the solar energy produced [cents/kWh]
							Can be a numpy array of length L, in case a dynamic price is applied.

		L                   The project term in years
		ResaleValue         The resale value of the PV system at the end of L years.

	RETURNS:
		cashflows               Annual cash flows starting with year 0 as the investment year.
		AnnualProduction_kWh    Annual solar PV energy produced in kWh. Starts with year 1.
	"""
	# ....................................................................... 1) Calculate annual solar production [kWh]
	AnnualProduction_kWh = calc_annualprod(L, deg_early, deg_cons, Sys_kW, solar_yield, start_at_zero=False)

	# ....................................................................... 2) Setup cashflows
	# a) Investment costs & solar revenue
	cashflows = np.append(-InvCosts, AnnualProduction_kWh*energy_price/100)
	# todo - add maintenance costs to AnnualProd

	# b) resale / end of life
	if ResaleValue != 0.0:
		#cashflows = np.append(cashflows, ResaleValue)
		cashflows[-1] += ResaleValue

	return cashflows, AnnualProduction_kWh


def calc_annualprod(L, deg_early=0.0, deg_cons=0.8, sys_kW=1.0, solar_yield=1.0, start_at_zero=True):
	"""Calculates the annual production of your PV system, considering degradation.

	PARAMETERS:
		L               Period in years

		deg_early       (Optional; defaults to 0) Early degradation [%]. Decrease in output after the first year. The
						default behavior doesn't simulate early degradation, and uses constant degradation.

		deg_cons        (Optional and recommended; defaults to 0.8) Steady degradation rate [%], after early
						degradation, if any.

		sys_kW*         (Optional; defaults to 1.0) System size in kW DC.
		solar_yield*    (Optional; defaults to 1.0) Solar yield in kWh(AC) per kW(DC).

		*Leave both to default values to get the nominal capacity in per-unit.

		start_at_zero   (Optional; defaults to True). If True, then the beginning of an investment timeline is
						considered, with a value of 0.0. Else, the first value represents the production/per unit
						capacity at year 1.


	RETURN:
		Numpy array of the per-unit capacity        (default sys_kW, solar_yield)
		                   capacity                 (actual sys_kW)
		                   annual production        (actual sys_kW and solar_yield)

	"""
	prod = np.ones(L)
	deg_cons /= 100
	deg_early /= 100

	if deg_early == 0.0:
		reduction = np.linspace(0, (L-1)*deg_cons, num=L)
		prod -= reduction

	else:
		reduction = np.linspace(0, (L-2)*deg_cons, num=(L-1))

		prod[1:] -= deg_early
		prod[1:] -= reduction

	results = sys_kW * solar_yield * prod

	if start_at_zero:
		return np.concatenate(([0], results))
	else:
		return results


def calc_resale(N, Sys_kW, solar_yield, energy_price, discount=0.1, useful_life=25, deg_early=0.0, deg_cons=0.8):
	"""Calculates the resale value of the PV system after N years, given useful life and factors to estimate the
	value of the generated solar energy.

	Method:
		The resale value of the PV system is estimated by discounting the remaining revenue in its useful life to the
		year of resale.

	PARAMETERS:
		N               Age of PV system at resale.
		useful_life     (Optional; default to 25 years) Years of useful life/warranty of PV system.
		discount        Discount rate to valuate the remaining future earnings.

		Sys_kW          The capacity of the PV system in kW DC.
		solar_yield     Solar yield in annual kWh(AC)/kW(DC)
		energy_price    The price of the solar-generated electricity [cents/kWh]

		deg_early       (Optional; defaults to 0) Argument to calc_annualprod()
		deg_cons        (Optional; defaults to 0.8) Argument to calc_annualprod()

	RETURNS:
		Resale value in dollars.
	"""
	remaining_yrs = useful_life-N
	remaining_cashflows  = calc_annualprod(useful_life, deg_early, deg_cons, Sys_kW, solar_yield)[N+1:] * \
	                       energy_price/100

	resale_value = round(npf.npv(discount, remaining_cashflows), 0)

	return resale_value


def plot_Price_vs_irr(PPA_prices, irrs, mark_tariffs=None, save_as=False, show_plot=True, **kwargs):
	"""Plot of PPA price vs IRR"""
	def_kws = {
		'title': 'Internal Rate of Return vs. Solar Energy Price',
		'title_kw': {'fontsize': 15},
		'ylabel': 'IRR (%)',
		'ylabel_kw': {'fontsize': 13},
		'yticks_kw': {'fontsize': 12},
		'xlabel': 'PPA energy price (US¢/kWh)',
		'xlabel_kw': {'fontsize': 13},
		'xticks_kw': {'fontsize': 12},
		'figsize': (8,6),
		'grid': True,
	}
	kwargs.update({key: val for key, val in def_kws.items() if key not in kwargs})

	# .................................................................................... a) Plot IRR vs PPA price
	plt.figure(figsize=kwargs['figsize'])
	plt.plot(PPA_prices, irrs * 100, color='#E67E22')

	# .................................................................................... b) Plot tariffs
	if mark_tariffs:
		tariff_y = 12.5

		for idx, keys in enumerate(mark_tariffs.keys()):
			tariff_lb, tariff_ub = mark_tariffs[keys]

			plt.plot(mark_tariffs[keys], np.array([tariff_y, tariff_y]), ls='--', lw=1, color='k')
			plt.plot(mark_tariffs[keys], np.array([tariff_y, tariff_y]), 'ko', ms=4)

			plt.text(0.2 * (tariff_ub - tariff_lb) + tariff_lb, tariff_y - 0.9, '{} tariffs'.format(keys))
			tariff_y += 2.5

	# .................................................................................... c) polish

	ax = plt.gca()
	ax = basic_plot_polishing(ax, **kwargs)

	plot_exit(ax, save_as, show_plot)
	return


def plot_InvCost_sensitivity(Invcosts_USDperW, IRRs, pu_Invcosts, markitems=None, save_as=False, show_plot=True,
                             **kwargs):
	"""Plots sensitivity of IRR with respect to investment costs.

	The last elements of the vectors are assumed to be the non-subsidized values.

	markitems = (-5, 0, 10, 20) --> -5%, nominal, +10%, +20%
	"""
	def_kws = {
		'title'    : 'IRR Sensitivity to Investment Costs',
		'title_kw' : {'fontsize': 15},
		'ylabel'   : 'IRR (%)',
		'ylabel_kw': {'fontsize': 13},
		'yticks_kw': {'fontsize': 12},
		'xlabel'   : 'Specific investment costs (USD/Wdc)',
		'xlabel_kw': {'fontsize': 13},
		'xticks_kw': {'fontsize': 12},
		'figsize'  : (8, 6),
		'grid'     : True,
	}
	kwargs.update({key: val for key, val in def_kws.items() if key not in kwargs})
	# .................................................................................... a) plot IRR vs inv costs
	plt.figure(figsize=kwargs['figsize'])
	plt.plot(Invcosts_USDperW, IRRs * 100, color='#E67E22')


	# .................................................................................... b) marks
	if markitems:
		marked_pu = np.array(markitems) / 100 + 1
		# indices of marked items in x-axis vec
		idxs = np.fromiter((np.where(pu_Invcosts == _pu)[0][0] for _pu in marked_pu), dtype='int')
		# Append unsubsidized
		idxs = np.append(idxs, -1)

		marked_Inv = np.fromiter((Invcosts_USDperW[idx] for idx in idxs), dtype='f8')
		marked_IRR = np.fromiter((IRRs[idx]*100 for idx in idxs), dtype='f8')

		plt.plot(marked_Inv, marked_IRR, 'o', color='#E67E22')

		# Add text
		for idx in range(marked_Inv.shape[0]):
			_x = marked_Inv[idx] + 0.02
			_y = marked_IRR[idx]

			try:
				mark = markitems[idx]
				if mark == 0:
					_text = 'nominal'
				elif mark > 0:
					_text = '+{}%'.format(mark)
				else:
					_text = '-{}%'.format(mark)
			except IndexError:
				_text = 'unsubsidized'
				_x -= 0.135
				_y -= 0.08

			plt.text(_x, _y, _text)

	# .................................................................................... c) polish
	ax = plt.gca()

	# set x-labels into 2-decimals
	locs, _ = plt.xticks()
	labels = ['{:.2f}'.format(itm) for itm in locs]
	plt.xticks(locs, labels)


	ax = basic_plot_polishing(ax, **kwargs)

	plot_exit(ax, save_as, show_plot)
	return


def plot_cashflows_solaronly(cashflows, resale_value, save_as=False, show_plot=True, **kwargs):
	def_kws = {
		'title'    : 'Revenue from generated solar',
		'title_kw' : {'fontsize': 15},
		'ylabel'   : 'kUSD',
		'ylabel_kw': {'fontsize': 13},
		'yticks_kw': {'fontsize': 12},
		'xlabel'   : 'Year',
		'xlabel_kw': {'fontsize': 13},
		'xticks_kw': {'fontsize': 12},
		'figsize'  : (8, 6),
		'grid': {'axis': 'y'},
		'ylims': (0, 55),
	}
	kwargs.update({key: val for key, val in def_kws.items() if key not in kwargs})
	# .................................................................................... a) plot solar revenue
	plt.figure(figsize=kwargs['figsize'])
	x_locs = list(range(len(cashflows)))

	# Get solar rev only
	solar_rev = np.copy(cashflows[1:])
	solar_rev[-1] -= resale_value
	solar_rev = solar_rev/1000

	plt.bar(x_locs[1:], solar_rev, color='#E67E22')

	# .................................................................................... b) polish
	ax = plt.gca()

	plt.xticks(x_locs, x_locs)

	ax = basic_plot_polishing(ax, **kwargs)
	plot_exit(ax, save_as, show_plot)
	return


def plot_cashflows_all(cashflows, resale_value, save_as=False, show_plot=True, **kwargs):
	def_kws = {
		'title'    : 'Cash Flows',
		'title_kw' : {'fontsize': 15},
		'ylabel'   : 'kUSD',
		'ylabel_kw': {'fontsize': 13},
		'yticks_kw': {'fontsize': 12},
		'xlabel'   : 'Year',
		'xlabel_kw': {'fontsize': 13},
		'xticks_kw': {'fontsize': 12},
		'figsize'  : (8, 6),
		'grid': {'axis': 'y'},
		'ylims': (-300, 300),
		'legend': True,
	}
	kwargs.update({key: val for key, val in def_kws.items() if key not in kwargs})
	# .................................................................................... a) plot
	plt.figure(figsize=kwargs['figsize'])
	x_locs = list(range(len(cashflows)))

	# -------------------------------------- i) solar rev only
	solar_rev = np.copy(cashflows[1:])
	solar_rev[-1] -= resale_value
	solar_rev = solar_rev/1000

	plt.bar(x_locs[1:], solar_rev, color='#E67E22', label='solar revenue')

	# -------------------------------------- ii) investment cost
	plt.bar(x_locs[0], cashflows[0]/1000, color='#B03A2E', label='investment costs')

	# -------------------------------------- iii) resale
	plt.bar(x_locs[-1], resale_value/1000, color='#2471A3', bottom=solar_rev[-1], label='system sale')


	# .................................................................................... b) polish
	ax = plt.gca()

	plt.xticks(x_locs, x_locs)

	ax = basic_plot_polishing(ax, **kwargs)
	plot_exit(ax, save_as, show_plot)
	return
































