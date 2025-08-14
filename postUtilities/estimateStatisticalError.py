import numpy as np
from scipy.optimize import minimize_scalar

def estimate_statistical_error(signal, dt):
    """
    Estimates the statistical error of a time signal based on the methods
    proposed by Mockett, Knacke, and Thiele (2010).

    The function calculates the empirical error trend by windowing the data,
    then fits an analytical white-noise error model to this trend to find
    an effective bandwidth 'B'. This 'B' is then used to estimate the
    statistical error for the mean and standard deviation of the entire signal.

    Args:
        signal (np.ndarray): A 1D numpy array representing the time signal.
        dt (float): The time step (delta t) between signal samples in seconds.

    Returns:
        dict: A dictionary containing the results:
              'total_mean': The mean of the entire signal.
              'total_std': The standard deviation of the entire signal.
              'optimal_B': The fitted effective bandwidth (Hz).
              'error_mean': The final estimated normalized error for the mean.
              'error_std': The final estimated normalized error for the stdev.
              'mean_95_confidence_interval': A tuple with the (lower, upper)
                                            bounds of the 95% confidence interval
                                            for the true mean.
              'empirical_data': A dict with 'window_times' and 'errors' for plotting.
    """
    if not isinstance(signal, np.ndarray):
        signal = np.array(signal)

    n_total = len(signal)
    if n_total < 10:
        raise ValueError("Signal is too short for a reliable statistical analysis.")

    t_total = n_total * dt
    mean_total = np.mean(signal)
    std_total = np.std(signal)

    if np.isclose(mean_total, 0):
        print("Warning: The mean of the signal is close to zero. "
              "Normalized error for the mean is not meaningful.")
        # We can still proceed to calculate B and error on std dev
        # but the error on the mean will be infinite.
        error_mean = np.inf
    
    # --- 1. Calculate Empirical Error Trend (Paper, Eq. 3) ---
    
    # Define a range of window sizes (Tw) to analyze.
    # We use a log scale to efficiently cover different time scales.
    # The paper suggests avoiding large Tw to have enough windows for stats.
    # We will go up to T/4 to ensure at least 4 windows.
    min_window_points = 10
    max_window_points = n_total // 4
    n_window_sizes = 30
    
    window_points = np.logspace(
        np.log10(min_window_points),
        np.log10(max_window_points),
        num=n_window_sizes
    ).astype(int)
    window_points = np.unique(window_points) # Remove duplicates

    empirical_errors = []
    valid_window_times = []

    for n_w in window_points:
        # Number of full windows that fit in the signal
        num_windows = n_total // n_w
        if num_windows < 4: # Need enough windows for a statistical sample
            continue\
            
        # Truncate signal to fit an integer number of windows
        truncated_signal = signal[:num_windows * n_w]
        
        # Reshape into a 2D array: (num_windows, window_size)
        windows = truncated_signal.reshape((num_windows, n_w))
        
        # Calculate the mean of each window
        window_means = np.mean(windows, axis=1)
        
        # The error is the standard deviation of the window means, normalized
        # by the total mean (as per the paper's definition of normalized error).
        std_of_means = np.std(window_means, ddof=1)
        
        if not np.isclose(mean_total, 0):
            error = std_of_means / abs(mean_total)
            empirical_errors.append(error)
            valid_window_times.append(n_w * dt)

    empirical_errors = np.array(empirical_errors)
    valid_window_times = np.array(valid_window_times)

    # --- 2. Fit Analytical Model to Find B (Paper, Eq. 4) ---
    
    # Analytical formula for error on the mean: ε ≈ 1 / sqrt(2 * B * T)
    def analytical_error_mean(B, T):
        if B <= 0:
            return np.inf
        return 1.0 / np.sqrt(2 * B * T)

    # Objective function to minimize: sum of squared log errors
    # Using log error gives better weighting to points across magnitudes
    def objective_function(B):
        log_analytical = np.log(analytical_error_mean(B, valid_window_times))
        log_empirical = np.log(empirical_errors)
        return np.sum((log_analytical - log_empirical)**2)

    # Find the optimal B that minimizes the objective function.
    # The search for B is bounded between a small number and the Nyquist frequency.
    nyquist_freq = 1.0 / (2 * dt)
    # Using minimize_scalar as it's efficient for one variable
    result = minimize_scalar(
        objective_function,
        bounds=(1e-9, nyquist_freq),
        method='bounded'
    )
    optimal_B = result.x

    # --- 3. Calculate Final Error and Confidence Interval ---
    
    # Use the optimal B to estimate the error for the *entire* signal length
    if not np.isclose(mean_total, 0):
        final_error_mean = analytical_error_mean(optimal_B, t_total)
    else:
        final_error_mean = np.inf
        
    # The paper uses the same B to estimate the error on the standard deviation
    # Analytical formula for error on the stdev: ε ≈ 1 / sqrt(4 * B * T)
    final_error_std = 1.0 / np.sqrt(4 * optimal_B * t_total)

    # Calculate the 95% confidence interval for the mean (Paper, Eq. 8)
    # True Mean is likely in [μ * (1 - 2ε), μ * (1 + 2ε)]
    lower_bound = mean_total * (1 - 2 * final_error_mean)
    upper_bound = mean_total * (1 + 2 * final_error_mean)

    return {
        'total_mean': mean_total,
        'total_std': std_total,
        'optimal_B': optimal_B,
        'error_mean': final_error_mean,
        'error_std': final_error_std,
        'mean_95_confidence_interval': (lower_bound, upper_bound),
        'empirical_data': {
            'window_times': valid_window_times,
            'errors': empirical_errors
        }
    }
def calculateAcrossTime(time, data,dt):
    ciUpperArray = []
    ciLowerArray = []
    avgDataArray = []
    for i in range(len(data)):
        try:
            results = estimate_statistical_error(data[0:i],dt=dt)
        except:
            continue
        ci = results['mean_95_confidence_interval']
        avg_data = results['total_mean']
        ciUpperArray.append(ci[1])
        ciLowerArray.append(ci[0])
        avgDataArray.append(avg_data)
    statStartTime = time[-1]-(len(ciUpperArray)*dt)
    statTime = np.linspace(statStartTime,time[-1],len(ciUpperArray))
    return statTime, avgDataArray,ciUpperArray,ciLowerArray


# if __name__ == '__main__':
#     # --- Create a sample signal for demonstration ---
#     # This signal is similar to the real-world data in the paper:
#     # a dominant periodic component (sine wave) with broadband random noise.
    
#     dt = 0.001  # Time step in seconds
#     T = 100     # Total time in seconds
#     t = np.arange(0, T, dt)
    
#     # Signal components
#     mean_val = 1.0
#     noise_std = 0.5
#     sine_amp = 0.25
#     sine_freq = 2.0 # Hz
    
#     # Generate the signal
#     noise = np.random.normal(0, noise_std, size=t.shape)
#     periodic = sine_amp * np.sin(2 * np.pi * sine_freq * t)
#     sample_signal = mean_val + noise + periodic
    
#     print(f"Generated a sample signal with {len(sample_signal)} points.")
#     print("-" * 30)

#     # --- Run the analysis ---
#     try:
#         results = estimate_statistical_error(sample_signal, dt)

#         # --- Print the results ---
#         print("Statistical Error Estimation Results:")
#         print(f"  - Signal Mean: {results['total_mean']:.4f}")
#         print(f"  - Signal Std Dev: {results['total_std']:.4f}")
#         print(f"  - Fitted Bandwidth (B): {results['optimal_B']:.2f} Hz")
#         print("\n--- Final Error Estimates ---")
#         print(f"  - Normalized Error on Mean (ε_μ): {results['error_mean']:.4f} "
#               f"({results['error_mean']*100:.2f}%)")
#         print(f"  - Normalized Error on Std Dev (ε_σ): {results['error_std']:.4f} "
#               f"({results['error_std']*100:.2f}%)")
        
#         ci = results['mean_95_confidence_interval']
#         print(f"\n  -> 95% Confidence Interval for the True Mean: "
#               f"[{ci[0]:.4f}, {ci[1]:.4f}]")

#         # --- Plot the results for visualization (highly recommended) ---
#         # This plot replicates the style of Fig. 1 and Fig. 4 in the paper
#         plt.figure(figsize=(10, 6))
        
#         # Plot the empirical data points
#         emp_data = results['empirical_data']
#         plt.loglog(
#             emp_data['window_times'], emp_data['errors'],
#             'x', label='Empirical Error (Eq. 3)', color='orangered'
#         )
        
#         # Plot the fitted analytical curve
#         T_fit = np.logspace(np.log10(emp_data['window_times'][0]), np.log10(T), 100)
#         B_fit = results['optimal_B']
#         error_fit = 1.0 / np.sqrt(2 * B_fit * T_fit)
#         plt.loglog(
#             T_fit, error_fit,
#             '-', label=f'Fitted Analytical Model (B={B_fit:.2f} Hz)', color='dodgerblue'
#         )
        
#         plt.title('Statistical Error Estimation: Empirical vs. Fitted Model')
#         plt.xlabel('Sample Length T_w (seconds)')
#         plt.ylabel('Normalized Error ε[μ]')
#         plt.grid(True, which='both', linestyle='--', linewidth=0.5)
#         plt.legend()
#         plt.show()

#     except ValueError as e:
#         print(f"Error: {e}")

