import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
from collections import Counter
import streamlit as st
from PIL import Image
import urllib.request
import io

# Padroniza a sequência e verifica se é DNA:


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



def eh_dna(seq):
    if set(seq).issubset({"A", "C", "G", "T", "*"}):
        return True
    else:
        return False


# Análise da composição de bases da sequência:

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


# Alinhamento com Needleman-Wunsch:

def cria_matriz_subs():

     matriz = {"col": ["A", "C", "G", "T"],
                 "A": [1, -2, -1, -2],
                 "C": [-2, 1, -2, -1],
                 "G": [-1, -2, 1, -2],
                 "T": [-2, -1, -2, 1]}

     return matriz


def calcula_score(base1, base2, matriz_subs):

    j = matriz_subs["col"].index(base1)

    for base in "ACGT":
        if base2 == base:
            score = matriz_subs[base2][j]
            return score


def valor_maximo(base1, base2, lado, cima, diagonal):

    if (base1 == base2) and (diagonal >= lado) and (diagonal >= cima):
        return diagonal
    elif (base1 != base2) and (diagonal >= lado) and (diagonal >= cima):
        return diagonal
    elif (lado > cima) and (lado > diagonal):
        return lado
    else:
        return cima


def cria_caminho(base1, base2, lado, cima, diagonal):

    if (base1 == base2) and (diagonal >= lado) and (diagonal >= cima):
        return "\\"
    elif (base1 != base2) and (diagonal >= lado) and (diagonal >= cima):
        return "\\"
    elif (lado > cima) and (lado > diagonal):
        return "-"
    else:
        return "|"


def lcs(seq1, seq2, matriz_subs):

    pontuacao = []
    caminho = []
    g = -3

    for i in range(0, len(seq1)):
        pontuacao.append([0] * len(seq2))
        caminho.append([""] * len(seq2))

    for i in range(0, len(seq1)):
        pontuacao[i][0] = g * i
        caminho[i][0] = "|"
    for j in range(0, len(seq2)):
        pontuacao[0][j] = g * j
        caminho[0][j] = "-"
 
    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):  

            base1 = seq1[i]
            base2 = seq2[j] 
            s = calcula_score(base1, base2, matriz_subs)
        
            lado = pontuacao[i][j-1]
            cima = pontuacao[i-1][j]
            diagonal = pontuacao[i-1][j-1]
    
            pontuacao[i][j] = valor_maximo(base1, base2, lado + g, cima + g, diagonal + s)
            caminho[i][j] = cria_caminho(base1, base2, lado + g, cima + g, diagonal + s)

    return caminho


def gera_alinhamento(seq1, seq2, matriz_caminho, matriz_subs):
    ali_seq1 = ""
    ali_seq2 = ""
    g = -3
    match = 0
    mismatch = 0
    gap = 0
    score_final = 0

    i = len(seq1)-1
    j = len(seq2)-1

    while (i != 0) or (j != 0):
        s = calcula_score(seq1[i], seq2[j], matriz_subs)

        if matriz_caminho[i][j] == "\\" and seq1[i] == seq2[j]:
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            match += 1
            score_final += s
            i -= 1
            j -= 1

        elif matriz_caminho[i][j] == "\\" and seq1[i] != seq2[j]:
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            mismatch += 1
            score_final += s
            i -= 1
            j -= 1    
    
        elif matriz_caminho[i][j] == "-":
            ali_seq1 = " - " + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            gap += 1
            score_final += g
            j -= 1
    
        elif matriz_caminho[i][j] == "|":
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = " - " + ali_seq2
            gap += 1
            score_final += g
            i -= 1

    return match, mismatch, gap, score_final, ali_seq1, ali_seq2


# função principal:

def main():

    url = "https://raw.githubusercontent.com/vanleiko/dna-streamlit/main/src/dna-logo2.jpg"
    with urllib.request.urlopen(url) as i:
	    byteImg = io.BytesIO(i.read())
	    imagem = Image.open(byteImg)

    default_input1 = """AGTTCGCACGGTTA"""
    default_input2 = """AGATTCGTACTGTA"""
    help_text = "Insira uma sequência de nucleotídeos ou uma sequência no formato FASTA"

    st.image(imagem, use_column_width=True)
    st.title("Alinhamento global de sequências de DNA")

    st.sidebar.markdown("""<p style='text-align: justify'>
    <b>🧬 Sobre</b><br>
    Web App que faz o alinhamento global entre duas sequências de DNA.<br><br> 
    <b>⚙️ Como funciona</b><br>
    Alinhamento baseado no algoritmo de Needleman-Wunsch (1970), o qual busca pelo alinhamento ótimo 
    entre duas sequências, ou seja, aquele com a maior pontuação final.<br><br> 
    <b>🔢 Pontuação utilizada</b><br> 
    <i>> Match:</i> +1, bases iguais alinhadas.<br>
    <i>> Mismatch:</i> -1, alinhamento entre purina-purina ou pirimidina-pirimidina.<br> 
    <i>> Mismatch:</i> -2, alinhamento entre purina-pirimidina ou pirimidina-purina.<br> 
    <i>> Gap:</i> -3, deleção ou inserção de base.<br>
    A matriz de substituição para scoring de match e mismatch foi baseada no modelo K2P (Kimura, 1980).<br><br>
    <b>📊 Composição das sequências</b><br>
    Esse Web App também analisa a quantidade e a porcentagem de Adenina, Citosina, Guanina e 
    Timina e o Conteúdo GC de cada sequência.</b><br><br> 
    <b>📚 Referências</b><br>
    Needleman, S.B. e Wunsch, C.D. 1970. A general method applicable to the search for similarities in 
    the amino acid sequences of two proteins. J. Mol. Bio., 48:443-453.<br>
    Kimura, M. 1980. A simple method for estimating evolutionary rates of base substitutions through comparative
    studies of nucleotide sequences. J. Mol. Evol., 16:111-120.
    </p>""", unsafe_allow_html=True)

    st.subheader("**Insira abaixo as suas sequências de DNA:**")       
    input_seq1 = st.text_area(label=">>> Sequência 1:", 
                        value=default_input1, height=100, help=help_text)
    input_seq2 = st.text_area(label=">>> Sequência 2:", 
                        value=default_input2, height=100, help=help_text)

    if input_seq1 and input_seq2:
        matriz_subs = cria_matriz_subs()
        seq1 = formata_sequencia(input_seq1)
        seq2 = formata_sequencia(input_seq2)
        
        if eh_dna(seq1) and eh_dna(seq2):

            matriz_caminho = lcs(seq1, seq2, matriz_subs)
            match, mismatch, gap, score_final, ali_seq1, ali_seq2 = gera_alinhamento(seq1, seq2, matriz_caminho, matriz_subs)
            
            df_quantidade = df_quantidade_nucleotideos(seq1, seq2)
            df_porcentagem = df_porcentagem_nucleotideos(seq1, seq2)
            df_gc = conteudo_gc(seq1, seq2)

            grafico_quantidade = gera_grafico(df_quantidade, "Quantidade") 
            grafico_porcentagem = gera_grafico(df_porcentagem, "Porcentagem") 
            grafico_gc = gera_grafico_gc(df_gc)
            
            st.subheader("**# 1. Alinhamento global:**")
            st.markdown(f"(1) {ali_seq1}<br>(2) {ali_seq2}", unsafe_allow_html=True)
            st.markdown(f"""<p><i>Matches: {match}<br>Mismatches: {mismatch}<br>Gaps: {gap}</i><br>
                        <b>Pontuação final = {score_final}""", unsafe_allow_html=True)   
                             
            st.subheader("**# 2. Análise das bases:**")            
            selecao = st.selectbox("Selecione o tipo de análise:", ["Quantidade", "Porcentagem"])
            resposta = st.radio("Mostrar gráficos?", ["Sim", "Não"], index=1) 

            if selecao == "Quantidade":
                st.dataframe(df_quantidade)
                if resposta == "Sim":
                    st.pyplot(grafico_quantidade)

            elif selecao =="Porcentagem":
                st.dataframe(df_porcentagem)
                if resposta == "Sim":
                    st.pyplot(grafico_porcentagem)
                    
            st.subheader("**# 3. Conteúdo GC:**")
            st.dataframe(df_gc)
            if resposta == "Sim":
                st.pyplot(grafico_gc)

        else:
            st.write("*As sequências inseridas não são DNA =(*")
            st.write("*Insira novamente suas sequências*")  

        
if __name__ == "__main__":
    main()