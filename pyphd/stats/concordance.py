# Functions for concordance analysis.
# Copyright (C) 2018  Lisa Perus
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# long with this program.  If not, see <https://www.gnu.org/licenses/>.

# System imports
from matplotlib import pyplot as plt

# Third Party imports
import numpy as np
from scipy import stats


def correlation_plot(xdata, ydata, xlim=[0, 1], ylim=[0, 1],
                     dataset_names=["Dataset 1", "Dataset 2"], out_png=None):
    """ Create a linear regression plot

    Parameters
    ----------
    xdata: numpy.ndarray
        First set of data for linear regression.
    ydata: numpy.ndarray
        Second set of data for linear regression.
    xlim: int array
        Limits for plot x axis.
    ylim: int array
        Limits for plot y axis.
    dataset_names: str arry
        Names of the datasets used as axis labels.
    out_png: str, default None
        path to output png, saved if not set to None.

    Returns
    -------
    out_png: str
        path to output png.
    """

    # Modify xdata for least square fit under the form [[x 1]] so that
    # ydata = a * xdata with a =[[slope intercept]]
    xdata_leastq = np.vstack([xdata, np.ones(len(xdata))]).T

    # Least squares fit of the data
    slope, intercept = np.linalg.lstsq(xdata_leastq, ydata)[0]
    slope = float(slope)
    intercept = float(intercept)

    # Calculate Pearson correlation coefficient and the 2-tailed p-value
    pearson_coeff, pval = stats.pearsonr(xdata, ydata)
    pearson_coeff = round(pearson_coeff, 2)

    # Set the axis limits for the plot
    plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]], fontsize='16', lw='2')

    # Plot a scatter plot
    plt.scatter(xdata, ydata, s=[1 for x in range(len(xdata))])

    # Plot the least squares fit line
    plt.plot(xdata, [slope*x + intercept for x in xdata], 'k', lw='2')

    # Add labels
    plt.xlabel(dataset_names[0], fontsize='16')
    plt.ylabel(dataset_names[1], fontsize='16')

    # Change the font size
    plt.tick_params(axis='both', which='major', labelsize='16', width='2')

    # Add pvalue to the plot
    if pval < 0.01:
            pval_str = '(P<0.01)'
    else:
            pval_str = '(P=' + str(round(pval, 3)) + ')'

    # Add a text box with correlation coefficient and p-value
    # TODO: Pass the text position as an argument
    plt.text(
        (xlim[1] + xlim[0])/2,
        (ylim[1] + ylim[0])/2, 'r=' + str(pearson_coeff) + pval_str,
        fontsize='16')

    # Set border thickness
    ax = plt.gca()
    ax.spines['top'].set_linewidth('2')
    ax.spines['left'].set_linewidth('2')
    ax.spines['right'].set_linewidth('2')
    ax.spines['bottom'].set_linewidth('2')

    if out_png is not None:
        plt.savefig(out_png)
    plt.close()

    return out_png, pearson_coeff, pval


def bland_altman_plot(xdata, ydata, xlim=[0, 1], ylim=[-1, 1],
                      xtitle="Average", ytitle="Difference", title=None,
                      out_png=None):
    """ Create a Bland Altman plot.

    Parameters
    ----------
    xdata: numpy.ndarray
        First set of data for the plot.
    ydata: numpy.ndarray
        Second set of data for the plot.
    xlim: int array
        Limits for plot x axis.
    ylim: int array
        Limits for plot y axis.
    xtitle: str
        Title for x axis.
    ytitle: str
        Title for y axis.
    title: str, default None
        Title for the plot

    Returns
    -------
    out_png: str
        path to output png.
    outliers_under: int array
        array of indexes where points difference is under the limit of
        agreement.
    outliers_over: int array
        array of indexes where points difference is over the limit of agreement
    """

    # Calculate the difference
    difference = ydata - xdata

    # Calculate the average values of the two datasets for each data point
    average = (xdata + ydata)/2

    # Calculate the mean of the differences
    mean_difference = np.mean(difference)

    # Calculate the standard deviation of the difference
    std_difference = np.std(difference)

    # Calculate the upper and lower limits of the agreement (95% confidence).
    upper_limit = mean_difference + 1.96*std_difference
    lower_limit = mean_difference - 1.96*std_difference

    # Set points over or under the upper and lower limit
    outliers_under = np.where(difference < lower_limit)
    outliers_over = np.where(difference > upper_limit)

    # Set axis limits for agreement
    plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]],
             fontsize='16', lw='2')

    # Plot
    plt.scatter(average, difference, s=[1 for x in range(len(xdata))])

    # Add the mean, upper and lower levels of agreement to the plot
    plt.axhline(y=mean_difference, lw='2', color='k', ls='dashed')
    plt.axhline(y=upper_limit, lw='2', color='k', ls='dashed')
    plt.axhline(y=lower_limit, lw='2', color='k', ls='dashed')

    # Change labeling and font size
    plt.xlabel(xtitle, fontsize='16')
    plt.ylabel(ytitle, fontsize='16')
    plt.tick_params(axis='both', which='major', labelsize='16', width='2')

    # Set plot title
    if title is not None:
        plt.title(title)

    # Set border thickness
    ax1 = plt.gca()
    ax1.spines['top'].set_linewidth('2')
    ax1.spines['left'].set_linewidth('2')
    ax1.spines['right'].set_linewidth('2')
    ax1.spines['bottom'].set_linewidth('2')

    # Save the plot
    if out_png is not None:
        plt.savefig(out_png)

    # End the function
    return out_png, outliers_under, outliers_over


def compute_lin_coeff(xdata, ydata):
    """ Computes Lin's concordance coefficient.
    (Lin L.I-K (1989) A concordance correlation coefficient to evaluate
     reproducibility. Biometrics 45:255-268)

    Parameters
    ----------
    xdata: numpy.ndarray
        First set of data.
    ydata: numpy.ndarray
        Second set of data.

    Returns
    -------
    coeff: float or None
        Lin's coefficient. Can be put to None if computation failed.
    pearson_coeff: float
        Pearson coefficient used to compute Lin coefficient.
    bias_correction_factor: float
        Bias correction factor used to compute Lin coefficient.
    lo: float
        Lower confidence interval value for Pearson coefficient.
    hi: float
        Higher confidence interval value for Pearson coefficient.
    """

    # Check for NaN values
    nan_values = np.where((np.isnan(xdata) | np.isnan(ydata)) == True)
    nan_values = nan_values[0]
    if len(nan_values) > 0:
        print("[Warning] : NaN values detected in x or y array data.")
        print("[Warning] : Discarding NaN values")
        condition = ((~ np.isnan(xdata)) & (~ np.isnan(ydata)))
        print(condition)
        xdata = xdata[condition]
        ydata = ydate[condition]
    if len(xdata) == 0:
        print("[Warning] : No values in xdata, aborting computation.")
        coeff = None
    elif len(ydata) == 0:
        print("[Warning] : No values in ydata, aborting computation.")
        coeff = None
    else:
        # Calculate numerator
        std_x = np.std(xdata)
        std_y = np.std(ydata)
        # pearson_coeff, pval = stats.pearsonr(xdata, ydata)
        pearson_coeff, pval, lo, hi = pearsonr_ci(xdata, ydata)

        mean_x = np.mean(xdata)
        mean_y = np.mean(ydata)
        numerator = 2 * std_x * std_y
        denominator = std_x**2 + std_y**2 + (mean_x - mean_y)**2
        bias_correction_factor = numerator / denominator
        coeff = bias_correction_factor * pearson_coeff

    return coeff, pearson_coeff, bias_correction_factor, lo, hi


def pearsonr_ci(xdata, ydata, alpha=0.05):
    """Calculate Pearson correlation along with the confidence levels.

    Parameters
    ----------
    xdata: numpy.ndarray
        First set of data.
    ydata: numpy.ndarray
        Second set of data.
    alpha : float
      Significance level for confidence levels computation. 0.05 by default.

    Returns
    -------
    r: float
     Pearson's correlation coefficient
    pval : float
      The corresponding p value
    lo, hi : float
      The lower and upper bound of confidence intervals
    """
    r, p = stats.pearsonr(xdata, ydata)
    r_z = np.arctanh(r)
    se = 1 / np.sqrt(xdata.size - 3)
    z = stats.norm.ppf(1 - alpha/2)
    lo_z, hi_z = r_z - z*se, r_z + z*se
    lo, hi = np.tanh((lo_z, hi_z))
    return r, p, lo, hi
