import pandas as pd
from typing import Union, Optional

def summarize_categorical(series: pd.Series) -> pd.DataFrame:
    """
    Summarize a categorical (or numeric treated as categorical) pandas Series.

    Parameters
    ----------
    series : pd.Series
        Input series to summarize.

    Returns
    -------
    pd.DataFrame
        DataFrame with counts and percentages sorted by frequency.
    """
    counts = series.value_counts(dropna=False)
    percentages = (series.value_counts(normalize=True, dropna=False) * 100).round(2)
    return pd.DataFrame({
        'Count': counts,
        'Percentage (%)': percentages
    }).sort_values("Count", ascending=False)


def categorize_numeric(
    series: pd.Series,
    high_label: str,
    high_value: float,
    mid_label: str,
    low_label: str,
    low_value: float
) -> pd.DataFrame:
    """
    Categorize numeric values into bins: high, middle, low.

    Parameters
    ----------
    series : pd.Series
        Numeric input series.
    high_label : str
        Label for values above `high_value`.
    high_value : float
        Threshold for high values.
    mid_label : str
        Label for values between `low_value` and `high_value`.
    low_label : str
        Label for values below `low_value`.
    low_value : float
        Threshold for low values.

    Returns
    -------
    pd.DataFrame
        Summary of counts and percentages per category.
    """
    total = len(series)
    bins = {
        high_label: (series > high_value).sum(),
        mid_label: ((series >= low_value) & (series <= high_value)).sum(),
        low_label: (series < low_value).sum()
    }
    df_bins = pd.DataFrame([
        {"Category": label, "Count": count, "Percentage (%)": round((count / total) * 100, 2)}
        for label, count in bins.items()
    ])
    return df_bins


def analyze_column(
    df: pd.DataFrame,
    column: Union[str, int],
    high_label: Optional[str] = None,
    high_value: Optional[float] = None,
    mid_label: Optional[str] = None,
    low_label: Optional[str] = None,
    low_value: Optional[float] = None,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Analyze a single column of a DataFrame: summary stats + frequency table.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    column : str or int
        Column name or index to analyze.
    high_label, high_value, mid_label, low_label, low_value : optional
        Parameters for numeric categorization.
    verbose : bool, default True
        If True, print summary information.

    Returns
    -------
    pd.DataFrame
        Summary table of the column (categorical counts or numeric categories).

    Example:
    --------
          ********************* ANALYSIS OF HBA **********************

      Unique values: 3
      Range: 1 – 3
      Mean: 1.89
      Median: 2.0


      Table:
      +----+------------+---------+------------------+
      |    |  Category  |   Count |   Percentage (%) |
      +====+============+=========+==================+
      |  0 |     a      |       2 |            10.53 |
      +----+------------+---------+------------------+
      |  1 |     b      |      13 |            68.42 |
      +----+------------+---------+------------------+
      |  2 |     c      |       4 |            21.05 |
      +----+------------+---------+------------------+

      ____________________________________________________________
    """
    # Select column
    if isinstance(column, int):
        series = df.iloc[:, column]
        column_name = df.columns[column]
    else:
        series = df[column]
        column_name = column

    if verbose:
        print('\n' + f"{' ANALYSIS OF ' + column_name + ' ':*^60}\n")

    if pd.api.types.is_numeric_dtype(series):
        summary_stats = {
            "Unique values": series.nunique(),
            "Range": f"{round(series.min(), 2)} – {round(series.max(), 2)}",
            "Mean": round(series.mean(), 2),
            "Median": round(series.median(), 2),
        }
        if verbose:
            for k, v in summary_stats.items():
                print(f"{k}: {v}")
            print()

        if None not in (high_label, high_value, mid_label, low_label, low_value):
            result = categorize_numeric(series, high_label, high_value, mid_label, low_label, low_value)
        else:
            if verbose:
                print("⚠️ No reference values provided for categorization.")
            result = summarize_categorical(series)

    else:
        result = summarize_categorical(series)
        if verbose:
            print("Type: Categorical")
            print(f"Unique entries: {len(result)}")
            top_val = result.index[0]
            low_val = result.index[-1]
            print(f"Most frequent: {top_val} ({result.loc[top_val, 'Percentage (%)']}%)")
            print(f"Least frequent: {low_val} ({result.loc[low_val, 'Percentage (%)']}%)")

    if verbose:
        print("\nTable:")
        print(result.to_markdown(tablefmt="grid", stralign="center"))
        print("\n" + "_" * 60)

    return result
