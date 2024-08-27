import pandas as pd
import numpy as np
import glob
import altair as alt

def theoretical_avgN(mu, temperature, box_area, lambda_square):
    return box_area * np.exp(mu / temperature) / lambda_square

def read_simulation_data(filename):
    timesteps = []
    num_particles = []
    pressures = []

    with open(filename, 'r') as file:
        for line in file:
            parts = line.split(',')
            timestep = int(parts[0].split(':')[1].strip())
            number_of_particles = int(parts[4].split(':')[1].strip())
            pressure = float(parts[3].split(':')[1].strip())
            timesteps.append(timestep)
            num_particles.append(number_of_particles)
            pressures.append(pressure)
    return timesteps, num_particles, pressures

if __name__ == "__main__":
    box_length_x = 20
    box_length_y = 20
    box_area = box_length_x * box_length_y
    lambda_square = 1.0  # Assuming lambda^2 = 1 for simplicity

    data_records = []

    # Read files and organize data by temperature
    for filename in glob.glob("data_*_*.txt"):
        filename_parts = filename.split('_')
        mu = float(filename_parts[1])
        temperature = float(filename_parts[2].replace('.txt', ''))

        # Filter out temperatures lower than 1
        if temperature < 1:
            continue

        _, num_particles, pressures = read_simulation_data(filename)
        sim_avgN = np.mean(num_particles)
        theory_avgN = theoretical_avgN(mu, temperature, box_area, lambda_square)
        avg_pressure = np.mean(pressures)
        stddevN = np.std(num_particles)/ np.size(num_particles)**0.5  # Calculate the standard deviation for error bars
        stddevP = np.std(pressures)/np.size(num_particles)**0.5  # Calculate the standard deviation for error bars

        data_records.append({
            'mu': mu,
            'temperature': temperature,
            'sim_avgN': sim_avgN,
            'theory_avgN': theory_avgN,
            'avg_pressure': avg_pressure,
            'stddevN': stddevN  # Add standard deviation to the record
        })

    # Convert to DataFrame and save to CSV
    df = pd.DataFrame(data_records)
    df.sort_values(by=['temperature', 'mu'], inplace=True)  # Ensure data is sorted
    df.to_csv('docs/simulation_data_ideal.csv', index=False)

    # Create Altair plot
    selection = alt.selection_single(fields=['temperature'], bind=alt.binding_range(min=1, max=10, step=1), name='Temperature')

    base = alt.Chart(df).transform_filter(selection)

    # Plot for Number of Particles with Error Bars
    particles_chart = base.mark_line(point=True).encode(
        x=alt.X('mu:Q', title='Chemical Potential (μ)'),
        y=alt.Y('sim_avgN:Q', title='Average Number of Particles (N)'),
        color=alt.value('blue'),
        tooltip=['mu', 'sim_avgN', 'theory_avgN', 'temperature']
    ).properties(
        title='Simulated vs Theoretical Average Number of Particles'
    ).encode(
        yError='stddevN:Q'  # Add error bars for the standard deviation
    )

    theory_particles_chart = base.mark_line(strokeDash=[5, 5], color='green').encode(
        x='mu:Q',
        y='theory_avgN:Q',
        tooltip=['mu', 'sim_avgN', 'theory_avgN', 'temperature']
    )

    # Plot for Pressure
    pressure_chart = base.mark_line(point=True).encode(
        x=alt.X('mu:Q', title='Chemical Potential (μ)'),
        y=alt.Y('avg_pressure:Q', title='Average Pressure'),
        color=alt.value('red'),
        tooltip=['mu', 'avg_pressure', 'temperature']
    ).properties(
        title='Average Pressure vs Chemical Potential'
    )

    # Combine plots
    combined_chart = alt.vconcat(
        particles_chart + theory_particles_chart,
        pressure_chart
    ).add_selection(
        selection
    ).properties(
        title='Ideal Gas in the GCMC Simulation:'
    )

    # Save the plot as an HTML file
    combined_chart.save('docs/N_P_plot_ideal.html')
