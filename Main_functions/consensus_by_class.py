import pandas as pd

def generate_dfs():
    """
    Genera y retorna los 4 DataFrames de referencia (Lab3, AI, Swiss, pkcsm).
    """
    #admet_lab_3_ = pd.concat([admet_lab_3_ref, admet_lab_3_top, admet_lab_3_control])
    admet_ai_ = pd.concat([admet_ai_ref, admet_ai_top, admet_ai_control])
    swiss_admet_ = pd.concat([swiss_admet_ref, swiss_admet_top, swiss_admet_control])
    pkcsm_ = pd.concat([pkcsm_ref, pkcsm_top, pkcsm_control])

    # Renombrar columnas específicas
    admet_ai_.rename(columns={'DILI': 'DILI_AI'}, inplace=True)

    return admet_ai_, swiss_admet_, pkcsm_ #admet_lab_3_


def calcular_puntaje_admet_por_rubro(dfs_list, **kwargs):
    """
    Calcula el puntaje ADMET para cada DataFrame de la lista.
    kwargs permite activar/desactivar rubros: FQ, MedChem, Abs, Dis, Met, Exc, Tox.
    """
    results_list = []

    for df, name in zip(dfs_list, [ "ADMET AI", "Swiss ADMET", "PKCSM"]):#"ADMET Lab 3",
        results = df.apply(
            calcular_puntaje_admet,
            axis=1,
            return_breakdown=True,
            **kwargs
        )

        # Separar resultados
        df['ADMET_Score'] = results.apply(lambda x: x[0])
        df['ADMET_Norm_Score'] = results.apply(lambda x: x[1])
        df['ADMET_Breakdown'] = results.apply(lambda x: x[2])
        df['ADMET_Rank'] = df['ADMET_Score'].rank(ascending=False, method='min')

        # Mostrar tabla resumida
        print(f"\n🧪 {name} Top:")
        final_df = df.loc[:, ['name', 'ADMET_Norm_Score', 'ADMET_Rank']]\
                     .sort_values(by='ADMET_Norm_Score', ascending=False)
        
        final_df_mod = final_df.applymap(replace_point_with_comma)
        #final_df_mod.rename(columns={'ADMET_Norm_Score': 'ADMET Score',
        #                             'ADMET_Rank': 'ADMET Rank',
        #                             'name': 'ID'
        #                             }, inplace=True)
        display(final_df_mod)

        results_list.append(final_df[['name', 'ADMET_Norm_Score']])

    return results_list


def consensus_score(dfs):
    """
    Genera un consenso de puntajes entre varios DataFrames.
    """
    dfs = [df.copy() for df in dfs]  # evitar modificar originales

    # Renombrar columnas dinámicamente
    for i, df in enumerate(dfs, start=1):
        df.rename(columns={"ADMET_Norm_Score": f"ADMET Score {i}"}, inplace=True)

    # Merge progresivo
    df_merged = dfs[0]
    for df in dfs[1:]:
        df_merged = pd.merge(df_merged, df, on="name", how="inner")

    # Detectar columnas de scores
    score_cols = [col for col in df_merged.columns if col.startswith("ADMET Score ")]

    # Calcular consenso
    df_merged["admet_score_consensus"] = df_merged[score_cols].mean(axis=1)

    # Ranking consenso
    df_merged["admet_rank_consensus"] = (
        df_merged["admet_score_consensus"]
        .rank(ascending=False, method="dense")
        .astype(int)
    )

    return df_merged.sort_values("admet_rank_consensus").reset_index(drop=True)


def replace_point_with_comma(x):
    """Convierte punto decimal a coma en floats."""
    return str(x).replace('.', ',') if isinstance(x, float) else x
