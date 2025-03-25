import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.ar_model import AutoReg

def ar_from_csv(fp, fn, fs=1000, order=10, fractional_scaling=13):
    """
    Reads a CSV file, extracts time-series data, fits an AR model, and returns scaled AR coefficients.

    Parameters:
    fp (str): File path.
    fn (str): File name.
    fs (int, optional): Sampling frequency in Hz (default is 1000 Hz).
    order (int, optional): Order of the AR model (default is 10).
    fractional_scaling (int, optional): Scaling factor for fixed-point conversion (default is 13).

    Returns:
    np.ndarray: Scaled AR coefficients as int32.
    """
    # Read the CSV file
    tbl = pd.read_csv(f"{fp}/{fn}")

    # Compute timestamps and create a time-based index
    t = tbl["dataTimestamp"] / fs  # Convert to seconds
    tbl["time"] = pd.to_timedelta(t, unit="s")  # Convert to timedelta format
    TT = tbl.set_index("time")  # Set time as index

    # Extract the second column (assumed to be the signal of interest)
    y = TT.iloc[:, 1]

    # Plot the data
    plt.figure()
    plt.plot(y, label="data")
    plt.legend()
    plt.show()

    # Fit an AR model using the Yule-Walker method
    ARmdl = AutoReg(y, lags=order, old_names=False).fit()

    # Extract and scale AR coefficients
    ar_coeffs = ARmdl.params[1:] 
    formatted_ar_coeffs = np.round(ar_coeffs * 2**fractional_scaling).astype(np.int32)

    return formatted_ar_coeffs

# Example usage:
# coeffs = compute_ar_coeffs("/path/to/folder", "filename.csv")
# print(coeffs)