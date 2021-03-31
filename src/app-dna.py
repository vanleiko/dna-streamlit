import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
from collections import Counter
import streamlit as st
from PIL import Image
import urllib.request
import io


@st.cache
def formata_sequencia(sequencia): 
    sequencia = sequencia.upper()
        
    if sequencia[0] == ">":
        sequencia = sequencia.splitlines()
        sequencia = sequencia[1:]
        sequencia = "".join(sequencia)
        sequencia = "*" + sequencia
        
    else:
        sequencia = sequencia.splitlines()
        sequencia = "".join(sequencia)
        sequencia = "*" + sequencia 
        
    return sequencia


@st.cache
def eh_dna(seq):
    if set(seq).issubset({"A", "C", "G", "T", "*"}):
        return True
    else:
        return False


def conta_nucleotideo(seq):
    
    mapa = Counter(seq)

    adenina = mapa["A"]
    citosina = mapa["C"]
    guanina = mapa["G"]
    timina = mapa["T"]
    total = len(seq)-1
    quantidade = [adenina, citosina, guanina, timina, total]

    porc_adenina = (adenina / total) * 100
    porc_citosina = (citosina / total) * 100
    porc_guanina = (guanina / total) * 100
    porc_timina = (timina / total) * 100
    porc_total = porc_adenina + porc_citosina + porc_guanina + porc_timina
    porcentagem = [porc_adenina, porc_citosina, porc_guanina, porc_timina, porc_total]

    return quantidade, porcentagem


def df_quantidade_nucleotideos(seq1, seq2):   

    df = pd.DataFrame(columns=["Feature", "Sequência 1", "Sequência 2"])
  
    quantidade1, _ = conta_nucleotideo(seq1)
    quantidade2, _ = conta_nucleotideo(seq2)

    df["Feature"] = ["Adenina", "Citosina", "Guanina", "Timina", "Total"]
    df["Sequência 1"] = quantidade1
    df["Sequência 2"] = quantidade2

    return df


def df_porcentagem_nucleotideos(seq1, seq2):   

    df = pd.DataFrame(columns=["Feature", "Sequência 1", "Sequência 2"])
  
    _, porcentagem1 = conta_nucleotideo(seq1)
    _, porcentagem2 = conta_nucleotideo(seq2)

    df["Feature"] = ["Adenina", "Citosina", "Guanina", "Timina", "Total"]
    df["Sequência 1"] = porcentagem1
    df["Sequência 2"] = porcentagem2

    return df


def gera_grafico(df, selecao="Quantidade"):

    df_melt = pd.melt(df[:4], id_vars="Feature", var_name="Sequência", value_name=selecao)

    fig, ax = plt.subplots(figsize=(6,3))
    ax = sns.barplot(x='Feature', y=selecao, hue='Sequência', data=df_melt)
    ax.set_title(f"{selecao} das bases", fontweight="bold", fontsize=14)
    ax.set_xlabel("Base", fontsize=12)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.legend(bbox_to_anchor=(1.35, 1))
  
    if selecao == "Porcentagem":
        ax.set_ylim(0, 100)
        ax.set_ylabel("Porcentagem (%)", fontsize=12)    
    else:
        ax.set_ylabel("Quantidade", fontsize=12) 

    return fig


def conteudo_gc(seq1, seq2):
       
    df = pd.DataFrame(columns=["Feature", "Sequência", "Porcentagem"])

    _, porcentagem1 = conta_nucleotideo(seq1)
    _, porcentagem2 = conta_nucleotideo(seq2)

    gc1 = porcentagem1[2] + porcentagem1[1]
    gc2 = porcentagem2[2] + porcentagem2[1]

    df["Feature"] = ["GC", "GC"]
    df["Sequência"] = ["Sequência 1", "Sequência 2"]
    df["Porcentagem"] = [gc1, gc2]
  
    return df


def gera_grafico_gc(df):
 
    fig, ax = plt.subplots(figsize=(6,3))
    ax = sns.barplot(x="Feature", y="Porcentagem", hue="Sequência", data=df)
    ax.set_title(f"Conteúdo GC", fontweight="bold", fontsize=14)
    ax.set_ylabel("Porcentagem (%)", fontsize=12)
    ax.set_ylim(0, 100)
    ax.set_xlabel(None)
    plt.legend(bbox_to_anchor=(1.35, 1))

    return fig


def valor_maximo(c1, c2, lado, cima, diagonal):
    if (c1 == c2) and (diagonal+1 >= lado) and (diagonal+1 >= cima):
        return diagonal+1  
    elif (lado >= cima) and (lado >= diagonal+1):
        return lado  
    else:
        return cima


def acha_caminho(c1, c2, lado, cima, diagonal):
    if (c1 == c2) and (diagonal+1 >= lado) and (diagonal+1 >= cima):
        return "\\"  
    elif (lado >= cima) and (lado >= diagonal+1):
        return "-"
    else:
        return "|"


def lcs(seq1, seq2):
  pontuacao = []
  caminho = []

  # preencher a matriz
  for i in range(0, len(seq1)):
    pontuacao.append([0] * len(seq2))
    caminho.append([""] * len(seq2))

  # preencher a primeira linha com "-" e a primeira coluna com "|"
  for i in range(0, len(seq1)):
    caminho[i][0] = "|"
  for j in range(0, len(seq2)):
    caminho[0][j] = "-"

  # linha
  for i in range(1, len(seq1)):
    # coluna
    for j in range(1, len(seq2)):
      # devolver a valor maximo do L
      pontuacao[i][j] = valor_maximo(seq2[j], seq1[i], pontuacao[i][j-1], pontuacao[i-1][j], pontuacao[i-1][j-1])
      caminho[i][j] = acha_caminho(seq2[j], seq1[i], pontuacao[i][j-1], pontuacao[i-1][j], pontuacao[i-1][j-1])
  
  return pontuacao, caminho


def gera_alinhamento(seq1, seq2, pontuacao, caminho):
    ali_seq1 = ""
    ali_seq2 = ""

    # linhas
    i = len(seq1)-1
    # colunas
    j = len(seq2)-1

    while (i != 0) or (j != 0):
        if caminho[i][j] == "\\":
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            i -= 1
            j -= 1
    
        elif caminho[i][j] == "-":
            ali_seq1 = " - " + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            j -= 1
        
        else:
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = " - " + ali_seq2
            i -= 1

    matches = pontuacao[-1][-1]

    return matches, ali_seq1, ali_seq2  


def main():

    url = "https://raw.githubusercontent.com/vanleiko/dna-streamlit/main/src/dna-logo.jpg"
    with urllib.request.urlopen(url) as i:
	    byteImg = io.BytesIO(i.read())
	    imagem = Image.open(byteImg)

    default_input1 = """ATGGCAACGGGATCGTAAAGCAAGATTCCGACTATCGTAGCTAGCTTGGAAAA"""
    default_input2 = """TCAATCGATCGTAAAGCAGATTCCGACTAAAGTAGCTAGCTTGTAAAT"""
    help_text = "Insira uma sequência de nucleotídeos ou uma sequência no formato FASTA"

    st.image(imagem, use_column_width=True)
    st.title("Alinhamento global de sequências de DNA")

    st.markdown("""<p style='text-align: justify'>
    Web App que faz o alinhamento global de sequências de DNA, par a par, usando o algoritmo de 
    <b>Needleman-Wunsch</b>, um dos algoritmos mais utilizados para alinhamento de sequências biológicas.<br> 
    Foi considerado apenas <b>match</b> entre as bases das sequências, com pontuação +1, enquanto 
    <b>gap</b> foi pontuação zero. Mismatch não foi considerado.<br> 
    Esse Web App também analisa a quantidade e a porcentagem de <b>A C G T</b> e o <b>Conteúdo GC</b> 
    de cada sequência.<br></p>""", unsafe_allow_html=True)
    st.write("Obs: realizei algumas melhorias neste web app, o qual pode ser acessado [neste link.](https://share.streamlit.io/vanleiko/dna-streamlit/main/src/app-dna-v2.py)")
    
    
    

    st.subheader("**Insira abaixo as suas sequências de DNA:**")       
    seq1 = st.text_area(label=">>> Sequência 1:", 
                        value=default_input1, height=150, help=help_text)

    seq2 = st.text_area(label=">>> Sequência 2:", 
                        value=default_input2, height=150, help=help_text)

    if seq1 and seq2:

        seq1_formatada = formata_sequencia(seq1)
        seq2_formatada = formata_sequencia(seq2)
        
        if eh_dna(seq1_formatada) and eh_dna(seq2_formatada) :

            matriz_pontuacao, matriz_caminho = lcs(seq1_formatada, seq2_formatada)
            matches, ali_seq1, ali_seq2 = gera_alinhamento(seq1_formatada, seq2_formatada, matriz_pontuacao, matriz_caminho)
            
            df_quantidade = df_quantidade_nucleotideos(seq1_formatada, seq2_formatada)
            df_porcentagem = df_porcentagem_nucleotideos(seq1_formatada, seq2_formatada)

            grafico_quantidade = gera_grafico(df_quantidade, "Quantidade") 
            grafico_porcentagem = gera_grafico(df_porcentagem, "Porcentagem") 
            
            df_gc = conteudo_gc(seq1_formatada, seq2_formatada)
            grafico_gc = gera_grafico_gc(df_gc)
            
            st.subheader("**# 1. Alinhamento global:**")
            st.write("Total de matches:", matches)
            st.write("(1)", ali_seq1)
            st.write("(2)", ali_seq2)
         
            st.subheader("**# 2. Análise das bases:**")
            selecao = st.selectbox("Selecione o gráfico:", ["Quantidade", "Porcentagem"])
     
            if selecao == "Quantidade":
                st.dataframe(df_quantidade)
                st.pyplot(grafico_quantidade)

            elif selecao =="Porcentagem":
                st.dataframe(df_porcentagem)
                st.pyplot(grafico_porcentagem)
                  
            st.subheader("**# 3. Conteúdo GC:**")
            st.dataframe(df_gc)
            st.pyplot(grafico_gc)

        else:
            st.write("*As sequências inseridas não são DNA =(*")
            st.write("*Insira novamente suas sequências*")        


if __name__ == "__main__":
    main()

