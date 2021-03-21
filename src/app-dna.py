import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
from PIL import Image


imagem = Image.open("dna-logo.jpg")
st.image(imagem, use_column_width=True)


@st.cache
def formata_sequencia(sequencia): 
    sequencia = sequencia.upper()

    if sequencia[0] == ">":
        sequencia = sequencia.splitlines()
        sequencia = sequencia[1:]
        sequencia = "".join(sequencia)
    
    else:
        sequencia = sequencia.splitlines()
        sequencia = "".join(sequencia)
    
    return sequencia


def eh_dna(sequencia):
    if set(sequencia).issubset({"A", "C", "G", "T"}):
        return True
    else:
        return False


def conta_nucleotideos(sequencia, selecao):    
 
    total = len(sequencia)

    adenina = sequencia.count("A")
    citosina = sequencia.count("C")
    guanina = sequencia.count("G")
    timina = sequencia.count("T")
    total = len(sequencia)

    porc_adenina = (adenina / total) * 100
    porc_citosina = (citosina / total) * 100
    porc_guanina = (guanina / total) * 100
    porc_timina = (timina / total) * 100
    porc_total = porc_adenina + porc_citosina + porc_guanina + porc_timina

    df = pd.DataFrame({"Quantidade": [adenina, citosina, guanina, timina, total],
                        "Porcentagem (%)": [porc_adenina, porc_citosina, porc_guanina, porc_timina, porc_total]},
                        index=["Adenina", "Citosina", "Guanina", "Timina", "Total"])
      
    fig, ax = plt.subplots(figsize=(6,3))
    ax = sns.barplot(x=df.index[:4], y=selecao, data=df.iloc[:4])
    ax.set_title(f"{selecao} de nucleotídeos", fontweight="bold", fontsize=14)
    ax.set_ylabel(f"{selecao}", fontsize=12)
    if selecao == "Porcentagem (%)":
        ax.set_ylim(0, 100)
  
    return df, fig


def conteudo_gc(sequencia):
    
    gc = sequencia.count("G") + sequencia.count("C")
    porc_gc = (gc / len(sequencia)) * 100
    porc_at = 100 - porc_gc
     
    df = pd.DataFrame({"Porcentagem (%)": [porc_gc, porc_at]},
                        index=["GC", "AT"])
    
    fig, ax = plt.subplots(figsize=(6,3))
    ax = sns.barplot(x=df.index, y="Porcentagem (%)", data=df)
    ax.set_title(f"Conteúdo GC", fontweight="bold", fontsize=14)
    ax.set_ylabel("Porcentagem (%)", fontsize=12)
    ax.set_ylim(0, 100)
    
    return df, fig


def main():
    
    st.title("# Contando nucleotídeos no DNA")
    st.markdown("Web App que conta a quantidade e porcentagem de **A C G T** e o **Conteúdo GC** de uma sequência de **DNA**.")

    st.subheader("**Insira a sua sequência de DNA**")
    help_texto = "Insira uma sequência de nucleotídeos ou uma sequência no formato FASTA"
    sequencia = st.text_area(label="Insira abaixo a sua sequência e pressione Ctrl + ENTER", height=250, help=help_texto)

    if sequencia:
        sequencia_formatada = formata_sequencia(sequencia)
        
        if eh_dna(sequencia_formatada):
  
            st.subheader("**# 1. Análise dos nucleotídeos:**")
            selecao = st.selectbox("Selecione o gráfico:", ["Quantidade", "Porcentagem (%)"])
            df_nucleotideo, grafico_nucleotideo = conta_nucleotideos(sequencia_formatada, selecao) 
            st.dataframe(df_nucleotideo)
            st.pyplot(grafico_nucleotideo)

            st.subheader("**# 2. Conteúdo GC:**")
            df_gc, grafico_gc = conteudo_gc(sequencia_formatada) 
            st.dataframe(df_gc)
            st.pyplot(grafico_gc)

        else:
            st.write("*A sequência inserida não é um DNA =(*")
            st.write("*Insira uma sequência de DNA*")


if __name__ == "__main__":
    main()

