import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import streamlit as st
from PIL import Image


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


def calcula_porcentagem(nucleotideo, total):
    porcentagem = (nucleotideo / total) * 100
    
    return porcentagem


def conta_nucleotideos(sequencia, selecao="Quantidade"): 
    conta = Counter(sequencia) 
    total = len(sequencia)
 
    adenina = conta["A"]
    citosina = conta["C"]
    guanina = conta["G"]
    timina = conta["T"] 
 
    porc_adenina = calcula_porcentagem(adenina, total)
    porc_citosina = calcula_porcentagem(citosina, total)
    porc_guanina = calcula_porcentagem(guanina, total)
    porc_timina = calcula_porcentagem(timina, total)
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
    conta = Counter(sequencia) 
    total = len(sequencia)
    gc = conta["G"] + conta["C"]
    
    porc_gc = calcula_porcentagem(gc, total)
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

    imagem = Image.open("dna-logo.jpg")
    default_input = open("sars_cov_2.fasta", "r").read()
    help_text = "Insira uma sequência de nucleotídeos ou uma sequência no formato FASTA"

    st.image(imagem, use_column_width=True)

    st.title("Contando nucleotídeos no DNA")
    st.markdown("""
    Web App que conta a quantidade e porcentagem de **A C G T** e o **Conteúdo GC** de uma sequência de **DNA**.
    """)

    st.subheader("**Insira a sua sequência de DNA**")       
    sequencia = st.text_area(label="Insira abaixo a sua sequência e pressione Ctrl + ENTER", 
                            value=default_input, height=250, help=help_text)

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

