import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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
