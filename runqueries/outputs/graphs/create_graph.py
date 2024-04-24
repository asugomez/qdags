import pandas as pd
import matplotlib.pyplot as plt

def plot_from_csv(csv_file):
    # Read the CSV file using pandas
    df = pd.read_csv(csv_file, delimiter=';', header=None)
    
    # Get the title of the graph
    title = df.iloc[0, 0]
    print(title)
    
    # Select the column for the x-axis (first column of the second row)
    x_column = df.iloc[1:, 0]
    print(x_column)
    
    # Get the names of the curves (starting from the second column of the second row)
    curve_names = df.iloc[0, 1:]
    
    print(curve_names)
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.title(title)
    plt.xlabel(x_column.iloc[0])
    plt.ylabel("Time")

    # Plot the curves
    for name in curve_names:
        plt.plot(x_column, df[name][1:], label=name)

    plt.legend()
    plt.show()



# Ejemplo de uso
plot_from_csv("../partial/louds/backtracking/results.csv")
