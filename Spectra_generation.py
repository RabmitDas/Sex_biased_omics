import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool, cpu_count

# Function to generate Lorentzian peak
def lorentzian(x, x0, gamma, A):
    return A * gamma**2 / ((x - x0)**2 + gamma**2)

# Function to generate synthetic NMR spectrum
def generate_spectrum(chemical_shifts, intensities, intensity_variation=0.2):
    ppm_range = np.linspace(-1, 10, 10000)
    spectrum = np.zeros_like(ppm_range)
    for shift, intensity in zip(chemical_shifts, intensities):
        lorentzian_peak = lorentzian(ppm_range*3.14*500*2, shift*3.14*500*2, 0.1, intensity)
        spectrum += lorentzian_peak
    return ppm_range, spectrum

# Load metabolite information from CSV file
metabolite_info = pd.read_csv("C:/Users/rabmi/Urine Metabolome - Sheet1.csv")

output_dir = "C:/Users/rabmi"  # Directory to save output CSV files

spectra_data = []
# Create output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
output_file = "Indivisual_spectra.csv"
all_spectra_df = pd.DataFrame(columns=["ppm"])

for i in range(0,31):
    metabolite_name = metabolite_info.at[i, "Compounds"]
    
    # Define the Excel file path using the metabolite name
    peak_csv = pd.read_csv(f"C:/Users/rabmi/Compound NMR/{metabolite_name}.csv")
        
    # Generate synthetic data for concentrations, chemical shifts, and intensities of compounds
    chemical_shifts = peak_csv["ppm"].tolist()
    intensities = peak_csv["Intensity"].tolist()

    # Simulate NMR spectrum for the mixture
    x, spectrum = generate_spectrum(chemical_shifts, intensities)
    y = x/(2*3.14*500)
    spectra_data.append({"Compound": metabolite_name, "Spectrum": spectrum, "Chemical_Shift": x})
    # Plot the simulated NMR spectrum
    plt.figure(figsize=(10, 6))
    plt.plot(x, spectrum, color='blue', label='Simulated Spectrum')
    plt.xlabel('Chemical Shift')
    plt.ylabel('Intensity')
    plt.title('Simulated NMR Spectrum for '+ metabolite_name)
    plt.legend()
    plt.grid(True)
    plt.show() 

# Concatenate all spectra data into a DataFrame
all_spectra_df = pd.concat([pd.DataFrame(data["Spectrum"], columns=[data["Compound"]], index=data["Chemical_Shift"]) for data in spectra_data], axis=1)

transposed_df = all_spectra_df.transpose()
# Save spectrum data to CSV file
transposed_df.to_csv(output_file)

# Define a function to generate synthetic NMR data for a single mixture
def generate_mixture(args):
    i, num_mixtures, ppm_range, spectra_df, concentration_range_df = args
    concentrations_list = []
    mixture_spectrum = np.zeros(len(ppm_range))
    for j in range(1, len(spectra_df.columns)):
        compound_name = spectra_df.columns[j]
        min_concentration_str = concentration_range_df.loc[concentration_range_df['Compounds'] == compound_name, "Minimum Concentration (µM/mM creatinine)"].values[0]
        max_concentration_str = concentration_range_df.loc[concentration_range_df['Compounds'] == compound_name, "Maximum Concentration (µM/mM creatinine)"].values[0]
        min_concentration = float(min_concentration_str)
        max_concentration = float(max_concentration_str)
        concentration = np.random.uniform(min_concentration, max_concentration)
        concentrations_list.append(concentration)
        scaled_spectrum = spectra_df.iloc[:, j] * concentration
        mixture_spectrum += scaled_spectrum
    noise_level = 0.2
    mixture_spectrum += noise_level * np.random.normal(size=len(ppm_range))
    return i, pd.DataFrame({f'Sample_{i}': mixture_spectrum}), pd.DataFrame({f'Sample_{i}': concentrations_list})

# Define a function to save the DataFrame to CSV
def save_to_csv(data, filename):
    data.to_csv(filename, index=False)

if __name__ == '__main__':
    # Define parameters
    num_mixtures = 200000
    output_dir = "/home/group5/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read the CSV files
    spectra_df = pd.read_csv("/home/group5/myenv/CSV/Indivisual_spectra.csv")
    concentration_range_df = pd.read_csv("/home/group5/myenv/CSV/Urine Metabolome - Sheet1.csv")
    ppm_range = spectra_df.iloc[:, 0]

    # Create pool of worker processes
    pool = Pool(cpu_count())

    # Generate synthetic NMR data for each mixture using multiprocessing
    args = [(i, num_mixtures, ppm_range, spectra_df, concentration_range_df) for i in range(num_mixtures)]
    results = pool.map(generate_mixture, args)

    # Concatenate all spectra data into a DataFrame
    all_spectra_data = [result[1] for result in results]
    all_label_data = [result[2] for result in results]
    all_spectra_df = pd.concat(all_spectra_data, axis=1)
    all_label_df = pd.concat(all_label_data, axis=1)
    Data_df = all_spectra_df.transpose()
    Label_df = all_label_df.transpose()
    # Save the DataFrame to CSV files using multiprocessing
    pool.apply_async(save_to_csv, (Data_df, os.path.join(output_dir, 'Data.csv')))
    pool.apply_async(save_to_csv, (Label_df, os.path.join(output_dir, 'Label.csv')))

    # Close the pool and wait for all tasks to complete
    pool.close()
    pool.join()

    print("Job Done")
