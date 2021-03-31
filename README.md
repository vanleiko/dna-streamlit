# dna-streamlit

Web App que faz o alinhamento global par a par de sequências de DNA e analisa a quantidade e a porcentagem de bases e o conteúdo GC das sequências.

Para o alinhamento, foi utilizado o algoritmo de Needleman-Wunsch, que faz alinhamento ótimo entre 2 sequências.

Na [versão 1](https://share.streamlit.io/vanleiko/dna-streamlit/main/src/app-dna.py) do webapp, foi considerado apenas match (+1) e gap (0) para o cálculo da pontuação para o alinhamento entre as sequências.

Na [versão 2](https://share.streamlit.io/vanleiko/dna-streamlit/main/src/app-dna-v2.py), foi utilizada uma matriz de scoring para match e mismatch, segundo Kimura, 1980. Sendo que:

Match = +1

Mismatch do tipo transição = -1

Mismatch do tipo transversão = -2

Para gap, a penalidade foi de -3.
