import pandas as pd
import matplotlib.pyplot as plt

def plot_from_csv(csv_file, name_plot, name_x, file_output):
    # Leer datos desde el archivo CSV
    data = pd.read_csv(csv_file, delimiter=';')

    # Variables del archivo
    variables = data.columns[1:]  # Excluir la columna 'k'

    # Crear el gráfico
    plt.figure(figsize=(12, 7))

    for variable in variables:
        plt.plot(data['k'], data[variable], marker='o', label=variable)

    # Configuración del gráfico
    plt.xscale('log')  # Escala logarítmica en el eje x
    plt.xlabel(name_x)
    plt.ylabel('Time (ms)')
    plt.title(name_plot)
    plt.legend()
    plt.grid(True)


    # Colocar la leyenda al lado del gráfico
    plt.legend(bbox_to_anchor=(1.05, .75), loc='upper left', borderaxespad=0.)

    # Mostrar el gráfico
    plt.tight_layout(rect=[0, 0, 0.85, 1])

    # Guardar el gráfico
    plt.savefig(file_output)

    # Mostrar el gráfico
    plt.show()


for type_fun in [0]:#,1]:
    file = f"../partial/louds/backtracking/results-f{type_fun}.csv"
    plot_from_csv(file, "Gradual retrieval (LOUDS - backtracking)", "k results", f"graphPartialLoudsBacktracking{type_fun}.png")
    file = f"../partial/louds/nonFixedQueue/results-f{type_fun}.csv"
    plot_from_csv(file, "Gradual retrieval (LOUDS - nonFixedQueue)", "k results", f"graphPartialLoudsNonFixedQueue{type_fun}.png")
    file = f"../partial/dfuds/backtracking/results-f{type_fun}.csv"
    plot_from_csv(file, "Gradual retrieval (DFUDS - backtracking)", "k results", f"graphPartialDfudsBacktracking{type_fun}.png")
    # file = f"../partial/dfuds/nonFixedQueue/results-f{type_fun}.csv"
    # plot_from_csv(file, "Gradual retrieval (DFUDS - nonFixedQueue)", "k results", f"graphPartialDfudsNonFixedQueue{type_fun}.png")
    #
    # file = f"../ranked/louds/backtracking/results-f{type_fun}.csv"
    # plot_from_csv(file, "Ranked retrieval (LOUDS - backtracking)", "top k results", f"graphRankedLoudsBacktracking{type_fun}.png")
    # file = f"../ranked/louds/nonFixedQueue/results-f{type_fun}.csv"
    # plot_from_csv(file, "Ranked retrieval (LOUDS - nonFixedQueue)", "top k results",f"graphRankedLoudsNonFixedQueue{type_fun}.png")
    # file = f"../ranked/dfuds/backtracking/results-f{type_fun}.csv"
    # plot_from_csv(file, "Ranked retrieval (DFUDS - backtracking)", "top k results", f"graphRankedDfudsBacktracking{type_fun}.png")
    # file = f"../ranked/dfuds/nonFixedQueue/results-f{type_fun}.csv"
    # plot_from_csv(file, "Ranked retrieval (DFUDS - nonFixedQueue)", "top k results", f"graphRankedDfudsNonFixedQueue{type_fun}.png")
