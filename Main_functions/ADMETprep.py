from functools import reduce
import pandas as pd

def normalize_columns(df, column_groups):
    """
    Normaliza los nombres de columnas de un DataFrame.

    Reemplaza las variantes de nombres de propiedades por un nombre estándar,
    de acuerdo con el diccionario `column_groups`.

    Parámetros
    ----------
    df : pandas.DataFrame
        DataFrame con las columnas originales.
    column_groups : dict
        Diccionario donde las llaves son los nombres estandarizados
        y los valores son listas de nombres equivalentes.

    Retorna
    -------
    pandas.DataFrame
        DataFrame con las columnas renombradas.
    """
    rename_map = {}
    for standard_name, variants in column_groups.items():
        for col in variants:
            if col in df.columns:
                rename_map[col] = standard_name
    return df.rename(columns=rename_map)


def applied_normalize_columns(df, list_column_groups):
    """
    Aplica la normalización de columnas usando múltiples diccionarios.

    Útil cuando se tienen varios grupos de propiedades con diferentes
    conjuntos de equivalencias de nombres.

    Parámetros
    ----------
    df : pandas.DataFrame
        DataFrame con las columnas originales.
    list_column_groups : list[dict]
        Lista de diccionarios de mapeo de columnas.

    Retorna
    -------
    pandas.DataFrame
        DataFrame con las columnas renombradas de acuerdo a todos
        los diccionarios proporcionados.
    """
    df = df.copy()
    for column_groups in list_column_groups:
        df = normalize_columns(df, column_groups)
    return df


def add_db_column(df, db_name):
    """
    Agrega una columna 'DB' al DataFrame.

    Esto permite etiquetar el origen de los datos cuando se combinan
    múltiples bases de datos en un mismo análisis.

    Parámetros
    ----------
    df : pandas.DataFrame
        DataFrame de entrada.
    db_name : str
        Nombre de la base de datos a asignar.

    Retorna
    -------
    pandas.DataFrame
        DataFrame con la columna 'DB' añadida.
    """
    df = df.copy()
    df["DB"] = db_name
    return df


def split_by_property_groups(df, list_column_groups):
    """
    Divide un DataFrame en sub-DataFrames según grupos de propiedades ADMET.

    Cada sub-DataFrame contiene solo las columnas relevantes para un grupo
    de propiedades (ej. Absorción, Distribución, Metabolismo, etc.).

    Además, genera reportes sobre qué columnas fueron asignadas y cuáles
    quedaron sin agrupar.

    Parámetros
    ----------
    df : pandas.DataFrame
        DataFrame con los datos ADMET.
    list_column_groups : list[tuple]
        Lista de tuplas en la forma (nombre_grupo, dict_columnas),
        donde `dict_columnas` es un diccionario de equivalencias.

    Retorna
    -------
    dict[str, pandas.DataFrame]
        Diccionario de sub-DataFrames, indexados por el nombre del grupo.
    """
    df_normalized = df.copy()
    all_subdfs = {}

    # Normalizar columnas para todos los grupos
    for _, column_groups in list_column_groups:
        df_normalized = normalize_columns(df_normalized, column_groups)

    # Crear un sub-DataFrame por grupo
    total_pro = 0
    cols = []
    for group_name, column_groups in list_column_groups:
        columns_needed = ['name', 'orig_can_smiles', 'DB'] + list(column_groups.keys())
        columns_in_df = [c for c in columns_needed if c in df_normalized.columns]
        sub_df = df_normalized[columns_in_df].copy()
        sub_df['Class'] = group_name
        all_subdfs[group_name] = sub_df

        # Reporte de columnas asignadas
        print(f"--- {group_name} ({len(columns_in_df)}) ---")
        print("Columnas asignadas:", columns_in_df)
        cols.append(columns_in_df)
        total_pro += len(columns_in_df)-3
        print()

    print(f"Cantidad total de columnas: {len(df_normalized.columns)}, cantidad asignada {total_pro}")

    # Revisar columnas no asignadas
    todas_las_columnas = []
    for _, column_groups in list_column_groups:
        todas_las_columnas.extend(column_groups.keys())
    original_columns = df_normalized.columns.tolist()

    no_asignadas = [col for col in original_columns if col not in todas_las_columnas
                    and col not in ['name', 'orig_can_smiles', 'DB']]

    print("Columnas originales no asignadas a ningún grupo:")
    print(no_asignadas)

    return all_subdfs


def get_sub_dfs(df, list_column_groups, name):
    """
    Atajo para obtener los sub-DataFrames principales de propiedades ADMET.

    Separa un DataFrame normalizado en los grupos más comunes:
    FQ, MedChem, Absorption, Distribution, Metabolism, Excretion y Toxicity.

    Parámetros
    ----------
    df : pandas.DataFrame
        DataFrame con datos ADMET.
    list_column_groups : list[tuple]
        Lista de tuplas en la forma (nombre_grupo, dict_columnas).
    name : str
        Identificador del conjunto de datos (para el reporte).

    Retorna
    -------
    tuple
        Sub-DataFrames en el orden:
        (FQ, MedChem, Absorption, Distribution, Metabolism, Excretion, Toxicity)
    """
    sub_dfs = split_by_property_groups(df, list_column_groups)
    df_fq = sub_dfs['FQ']
    df_medchem = sub_dfs['MedChem']
    df_abs = sub_dfs['Absorption']
    df_dis = sub_dfs['Distribution']
    df_metab = sub_dfs['Metabolism']
    df_exc = sub_dfs['Excretion']
    df_tox = sub_dfs['Toxicity']

    print(f"====== {name} ======\n")
    return df_fq, df_medchem, df_abs, df_dis, df_metab, df_exc, df_tox


def display_min_max(df, list_of_properties):
    x = df[list_of_properties]
    a = x.min()
    b = x.max()
    return x, a, b


def assign_color_by_property(row, rules):
    '''
    elif rule["type"] == "categorical":
        yes, no = rule["thresholds"]
        if value == yes:
            return "🟢", "Verde"
        elif value == no:
            return "🔴", "Rojo"
        else:
            return "⚪", "Desconocido"
    '''
    prop = row["Property"]
    value = row["Value"]

    if prop not in rules:
        return "⚪", "Desconocido"

    rule = rules[prop]

    # -----------------------
    # Valores numéricos simples
    # -----------------------
    if rule["type"] == "numeric_up": #ascendente
        low, high = rule["thresholds"]
        if value <= low:
            return "🟢", f"Verde (óptimo: Valor <= {low})"
        elif value <= high:
            return "🟡", f"Amarillo (óptimo: Valor <= {low} )"
        else:
            return "🔴", f"Rojo (óptimo: Valor <= {low})"

    if rule["type"] == "numeric_down": #descendente
        high, low  = rule["thresholds"]
        if value >= high:
            return "🟢", f"Verde (Valor >= {high})"
        elif value >= low:
            return "🟡", f"Amarillo (Valor >= {high})"
        else:
            return "🔴", f"Rojo (Valor >= {high})"
    # -----------------------
    # Valores categóricos (flexibles)
    # -----------------------
    elif rule["type"] == "categorical":
        thresholds = rule["thresholds"]
        for categories, color in thresholds:
            if value in categories:
                if color == "Green":
                    return "🟢", f"Verde"
                elif color == "Yellow":
                    return "🟡", "Amarillo"
                elif color == "Red":
                    return "🔴", "Rojo"
        return "⚪", "Desconocido"

    # -----------------------
    # Rango definido como "óptimo"
    # -----------------------
    elif rule["type"] == "range":
        low, high = rule["green"]
        if low <= value <= high:
            return "🟢", f"Óptimo ({low}-{high})"
        elif value < low:
            return "🔴", f"Pequeño (óptimo: {low}-{high})"
        else:
            return "🔴", f"Grande (óptimo: {low}-{high})"

    return "⚪", "Desconocido"
