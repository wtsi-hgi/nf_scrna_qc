#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-06-29'
__version__ = '0.0.1'

import argparse
import warnings
import numpy as np
# from scipy import stats
import pandas as pd
import plotnine as plt9
import scanpy as sc


# NOTE: copied this code from scipy. It was introduced with 1.5.0, but due to a
# bug in scanpy dendrogram (1.5.1) scipy >= 1.5.0 is not compatible with
# scanpy.
# https://github.com/scipy/scipy/blob/v1.5.1/scipy/stats/stats.py#L2878-L3028
# START copied code -----------------------------------------------------------
def _contains_nan(a, nan_policy='propagate'):
    policies = ['propagate', 'raise', 'omit']
    if nan_policy not in policies:
        raise ValueError("nan_policy must be one of {%s}" %
                         ', '.join("'%s'" % s for s in policies))
    try:
        # Calling np.sum to avoid creating a huge array into memory
        # e.g. np.isnan(a).any()
        with np.errstate(invalid='ignore'):
            contains_nan = np.isnan(np.sum(a))
    except TypeError:
        # This can happen when attempting to sum things which are not
        # numbers (e.g. as in the function `mode`). Try an alternative method:
        try:
            contains_nan = np.nan in set(a.ravel())
        except TypeError:
            # Don't know what to do. Fall back to omitting nan values and
            # issue a warning.
            contains_nan = False
            nan_policy = 'omit'
            warnings.warn(
                "The input array could not be properly checked for nan "
                "values. nan values will be ignored.", RuntimeWarning
            )

    if contains_nan and nan_policy == 'raise':
        raise ValueError("The input contains nan values")

    return contains_nan, nan_policy


def _mad_1d(x, center, nan_policy):
    # Median absolute deviation for 1-d array x.
    # This is a helper function for `median_abs_deviation`; it assumes its
    # arguments have been validated already.  In particular,  x must be a
    # 1-d numpy array, center must be callable, and if nan_policy is not
    # 'propagate', it is assumed to be 'omit', because 'raise' is handled
    # in `median_abs_deviation`.
    # No warning is generated if x is empty or all nan.
    isnan = np.isnan(x)
    if isnan.any():
        if nan_policy == 'propagate':
            return np.nan
        x = x[~isnan]
    if x.size == 0:
        # MAD of an empty array is nan.
        return np.nan
    # Edge cases have been handled, so do the basic MAD calculation.
    med = center(x)
    mad = np.median(np.abs(x - med))
    return mad


def median_abs_deviation(x, axis=0, center=np.median, scale=1.0,
                         nan_policy='propagate'):
    r"""
    Compute the median absolute deviation of the data along the given axis.

    The median absolute deviation (MAD, [1]_) computes the median over the
    absolute deviations from the median. It is a measure of dispersion
    similar to the standard deviation but more robust to outliers [2]_.

    The MAD of an empty array is ``np.nan``.

    .. versionadded:: 1.5.0

    Parameters
    ----------
    x : array_like
        Input array or object that can be converted to an array.
    axis : int or None, optional
        Axis along which the range is computed. Default is 0. If None, compute
        the MAD over the entire array.
    center : callable, optional
        A function that will return the central value. The default is to use
        np.median. Any user defined function used will need to have the
        function signature ``func(arr, axis)``.
    scale : scalar or str, optional
        The numerical value of scale will be divided out of the final
        result. The default is 1.0. The string "normal" is also accepted,
        and results in `scale` being the inverse of the standard normal
        quantile function at 0.75, which is approximately 0.67449.
        Array-like scale is also allowed, as long as it broadcasts correctly
        to the output such that ``out / scale`` is a valid operation. The
        output dimensions depend on the input array, `x`, and the `axis`
        argument.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan.
        The following options are available (default is 'propagate'):

        * 'propagate': returns nan
        * 'raise': throws an error
        * 'omit': performs the calculations ignoring nan values

    Returns
    -------
    mad : scalar or ndarray
        If ``axis=None``, a scalar is returned. If the input contains
        integers or floats of smaller precision than ``np.float64``, then the
        output data-type is ``np.float64``. Otherwise, the output data-type is
        the same as that of the input.

    See Also
    --------
    numpy.std, numpy.var, numpy.median, scipy.stats.iqr, scipy.stats.tmean,
    scipy.stats.tstd, scipy.stats.tvar

    Notes
    -----
    The `center` argument only affects the calculation of the central value
    around which the MAD is calculated. That is, passing in ``center=np.mean``
    will calculate the MAD around the mean - it will not calculate the *mean*
    absolute deviation.

    The input array may contain `inf`, but if `center` returns `inf`, the
    corresponding MAD for that data will be `nan`.

    References
    ----------
    .. [1] "Median absolute deviation",
           https://en.wikipedia.org/wiki/Median_absolute_deviation
    .. [2] "Robust measures of scale",
           https://en.wikipedia.org/wiki/Robust_measures_of_scale

    Examples
    --------
    When comparing the behavior of `median_abs_deviation` with ``np.std``,
    the latter is affected when we change a single value of an array to have an
    outlier value while the MAD hardly changes:

    >>> from scipy import stats
    >>> x = stats.norm.rvs(size=100, scale=1, random_state=123456)
    >>> x.std()
    0.9973906394005013
    >>> stats.median_abs_deviation(x)
    0.82832610097857
    >>> x[0] = 345.6
    >>> x.std()
    34.42304872314415
    >>> stats.median_abs_deviation(x)
    0.8323442311590675

    Axis handling example:

    >>> x = np.array([[10, 7, 4], [3, 2, 1]])
    >>> x
    array([[10,  7,  4],
           [ 3,  2,  1]])
    >>> stats.median_abs_deviation(x)
    array([3.5, 2.5, 1.5])
    >>> stats.median_abs_deviation(x, axis=None)
    2.0

    Scale normal example:

    >>> x = stats.norm.rvs(size=1000000, scale=2, random_state=123456)
    >>> stats.median_abs_deviation(x)
    1.3487398527041636
    >>> stats.median_abs_deviation(x, scale='normal')
    1.9996446978061115

    """
    if not callable(center):
        raise TypeError("The argument 'center' must be callable. The given "
                        f"value {repr(center)} is not callable.")

    # An error may be raised here, so fail-fast, before doing lengthy
    # computations, even though `scale` is not used until later
    if isinstance(scale, str):
        if scale.lower() == 'normal':
            scale = 0.6744897501960817  # special.ndtri(0.75)
        else:
            raise ValueError(f"{scale} is not a valid scale value.")

    x = np.asarray(x)

    # Consistent with `np.var` and `np.std`.
    if not x.size:
        if axis is None:
            return np.nan
        nan_shape = tuple(item for i, item in enumerate(x.shape) if i != axis)
        if nan_shape == ():
            # Return nan, not array(nan)
            return np.nan
        return np.full(nan_shape, np.nan)

    contains_nan, nan_policy = _contains_nan(x, nan_policy)

    if contains_nan:
        if axis is None:
            mad = _mad_1d(x.ravel(), center, nan_policy)
        else:
            mad = np.apply_along_axis(_mad_1d, axis, x, center, nan_policy)
    else:
        if axis is None:
            med = center(x, axis=None)
            mad = np.median(np.abs(x - med))
        else:
            # Wrap the call to center() in expand_dims() so it acts like
            # keepdims=True was used.
            med = np.expand_dims(center(x, axis=axis), axis)
            mad = np.median(np.abs(x - med), axis=axis)

    return mad / scale
# END copied code -----------------------------------------------------------


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Calculate one, two or three MADs above median for qc metric \
            specified in key. Saves tsv file with calculated MADs and \
            distribution plot of qc key metric.
            """
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-qck', '--qc_key',
        action='store',
        dest='qc_key',
        default='pct_counts_gene_group__mito_transcript',
        help='QC key for MADs calculation.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output png file. Will have -mads-<qc_key> appended\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    output_file = options.of
    if output_file == '':
        output_file = 'mads-{}'.format(
            options.qc_key
        )

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)

    qc_keys = options.qc_key.split(',')
    df_dict = {
        'variable': [],
        'experiment_id': [],
        'cutoff_type': [],
        'cutoff': [],
        'cutoff_round': []
    }
    for qc_key in qc_keys:
        med = np.median(adata.obs[qc_key])
        # mad = stats.median_abs_deviation(adata.obs[qc_key])
        mad = median_abs_deviation(adata.obs[qc_key])
        for i in [1, 2, 3]:
            df_dict['variable'].append(qc_key)
            df_dict['experiment_id'].append('all_experiment_ids')
            df_dict['cutoff_type'].append('median+({}*mad)'.format(i))
            df_dict['cutoff'].append(med+(i*mad))
            df_dict['cutoff_round'].append(round(med+(i*mad)))
        for i in [1, 2, 3]:
            df_dict['variable'].append(qc_key)
            df_dict['experiment_id'].append('all_experiment_ids')
            df_dict['cutoff_type'].append('median-({}*mad)'.format(i))
            df_dict['cutoff'].append(med-(i*mad))
            df_dict['cutoff_round'].append(round(med-(i*mad)))

    df = pd.DataFrame(df_dict)
    df.to_csv(
        '{}.tsv'.format(output_file),
        sep='\t',
        index=False
    )

    for qc_key in qc_keys:
        gplt = plt9.ggplot(adata.obs, plt9.aes(x=qc_key))
        gplt = gplt + plt9.geom_histogram()
        filt = (df['variable'] == qc_key) & np.in1d(df['cutoff_type'].values, [
            'median+({}*mad)'.format(i) for i in [1, 2, 3]
        ])
        for cut in df.loc[filt, 'cutoff']:
            gplt = gplt + plt9.geom_vline(
                xintercept=cut,
                linetype="dashed",
                color="red"
            )
        gplt.save(
            '{}-{}.png'.format(output_file, qc_key),
            #dpi=300,
            width=5,
            height=4
        )


if __name__ == '__main__':
    main()
