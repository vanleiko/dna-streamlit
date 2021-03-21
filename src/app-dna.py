import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st


def formata_sequencia(sequencia):    
    if sequencia[0] == ">":
        sequencia = sequencia.splitlines()
        sequencia = sequencia[1:]
        sequencia = "".join(sequencia)
    
    else:
        sequencia = sequencia.splitlines()
        sequencia = "".join(sequencia)
    
    return sequencia


def conta_nucleotidos(sequencia):    
    sequencia = sequencia.upper()

    adenina = sequencia.count("A")
    citosina = sequencia.count("C")
    guanina = sequencia.count("G")
    timina = sequencia.count("T")

    df = pd.DataFrame({"Quantidade": [adenina, citosina, guanina, timina]},
                        index=["Adenina", "Citosina", "Guanina", "Timina"])

    fig, ax = plt.subplots(figsize=(6,3))
    ax = sns.barplot(x=df.index, y="Quantidade", data=df)
    ax.set_title(f"Contagem de nucleotídeos", fontweight="bold", fontsize=14)
    ax.set_ylabel("Quantidade", fontsize=12)
  
    return df, fig


def conteudo_gc(sequencia):
    sequencia = sequencia.upper()

    gc = sequencia.count("G") + sequencia.count("C")
    porcentagem_gc = (gc / len(sequencia)) * 100
    porcentagem_at = 100 - porcentagem_gc
     
    df = pd.DataFrame({"Porcentagem (%)": [porcentagem_gc, porcentagem_at]},
                        index=["GC", "AT"])
    
    fig, ax = plt.subplots(figsize=(6,3))
    ax = sns.barplot(x=df.index, y="Porcentagem (%)", data=df)
    ax.set_title(f"Conteúdo GC", fontweight="bold", fontsize=14)
    ax.set_ylabel("Porcentagem (%)", fontsize=12)
    ax.set_ylim(0, 100)
    
    return df, fig


def main():
    
    st.title("Contagem de Nucleotídeos e Conteúdo GC do DNA")
    st.markdown("Web App que conta a quantidade de **A C G T** e o **Conteúdo GC** de uma sequência de **DNA**.")

    st.subheader("**Insira sua sequência de DNA:**")
    seq_input = """>Meu_DNA
                    GAACACGTGGAGGCAAACAGGAAGGTGAAGAAGAACTTATCCTATCAGGACGGAAGGTGCTCGG
                    ATCTTCCTCGCGACTCTAAATTGCCCCCTCTGAGGTCAAGGAACACAAGATGGTTTTGGAAATG
                    TGAACCCATTATAACATCACCAGCATCGTGCCTGAAGCCATGCCTGCTGCCACCATGCCAGTCC"""
    sequencia = st.text_area("Insira uma sequência de nucleótideos ou formato FASTA", seq_input, height=300)
    sequencia_formatada = formata_sequencia(sequencia)
    st.write("Sequência inserida:", sequencia_formatada)
       
    st.subheader("**1. Quantidade de nucleotídeos:**")
    df_nucleotideo, grafico_nucleotideo = conta_nucleotidos(sequencia_formatada) 
    st.dataframe(df_nucleotideo)
    st.pyplot(grafico_nucleotideo)

    st.subheader("**2. Conteúdo GC:**")
    df_gc, grafico_gc = conteudo_gc(sequencia_formatada) 
    st.dataframe(df_gc)
    st.pyplot(grafico_gc)


if __name__ == "__main__":
    main()



