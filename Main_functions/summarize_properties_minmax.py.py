import pandas as pd
def summarize_property_ranges(df, id_column, props, group_dict, predictor, decimals=2):
    """
    Genera una tabla comparativa de propiedades por grupo de moléculas.

    Args:
        df (pd.DataFrame): DataFrame con propiedades.
        id_column (str): Columna que contiene los identificadores de moléculas.
        props (list): Lista de propiedades (columnas) a comparar.
        group_dict (dict): Diccionario con nombre del grupo y lista de IDs.
        decimals (int): Número de decimales para propiedades numéricas.

    Returns:
        pd.DataFrame: Tabla resumen con rangos o frecuencias por grupo.

    Ejemplo:
          Propiedad	Predictor	Total (60)	  Top 10	      Referencias	Control Positivo
        0	f_Csp3	  SwissADME	0.75 — 1.0	  0.75 — 0.94	  1.0	               0.9
    """
    summary_rows = []

    for prop in props:
        row = {'Propiedad': prop}

        for group_name, id_list in group_dict.items():
            subset = df[df[id_column].isin(id_list)]

            if subset.empty or prop not in subset.columns:
                row[group_name] = "N/A"
                continue

            values = subset[prop].dropna()

            if values.empty:
                row[group_name] = "N/A"
                continue

            # Datos numéricos
            if pd.api.types.is_numeric_dtype(values):
                min_val = round(values.min(), decimals)
                max_val = round(values.max(), decimals)
                if min_val == max_val:
                    row[group_name] = f"{min_val}"
                else:
                    row[group_name] = f"{min_val} — {max_val}"

            # Datos categóricos / booleanos
            else:
                counts = values.value_counts()
                counts_str = ', '.join([f"{str(k)} ({v})" for k, v in counts.items()])
                row[group_name] = counts_str

        summary_rows.append(row)
    table = pd.DataFrame(summary_rows)
    #insertar columna en df en posición 1
    table.insert(loc=1, column='Predictor', value=predictor)

    return table

def get_groups(general, top, ref, ctrl):
  group_dict = {
    'Total (60)': general['name'].tolist(),
    'Top 10': top['name'].tolist(),
    'Referencias': ref['name'].tolist(),
    'Control Positivo': ctrl['name'].tolist()
  }
  return group_dict


def applied_summary(df, df_top, df_ref, df_ctrl,
                    predictor, id_column='name'):
  df_group = get_groups (df, df_top, df_ref, df_ctrl)
  tabla = summarize_property_ranges(df, id_column, props=df.columns.drop(id_column).tolist(), group_dict=df_group, predictor=predictor)
  return tabla
