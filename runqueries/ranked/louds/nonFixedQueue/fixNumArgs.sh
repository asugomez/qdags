#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Uso: $0 <archivo>"
    exit 1
fi

archivo="$1"

# Verificar si el archivo existe
if [ ! -f "$archivo" ]; then
    echo "El archivo $archivo no existe."
    exit 1
fi

# Crear un archivo temporal para almacenar la salida
archivo_temporal=$(mktemp)

# Leer el archivo línea por línea y procesarlo
while IFS= read -r linea; do
    # Eliminar "./j4" del inicio de la línea y dividir los elementos por espacios
    elementos=($(echo "$linea" | sed 's/\.\// /' | awk '{for (i=2; i<=NF; i++) print $i}'))

    # Imprimir la línea original junto con los elementos al final
    printf "%s" "$linea" >> "$archivo_temporal"
    for ((i = 0; i < ${#elementos[@]}; i++)); do
        printf " %s" "${elementos[$i]}" >> "$archivo_temporal"
    done
    printf "\n" >> "$archivo_temporal"
done < "$archivo"

# Reemplazar el archivo original con el archivo temporal
mv "$archivo_temporal" "$archivo"

echo "Se han agregado los elementos al final de cada línea en el archivo $archivo."
