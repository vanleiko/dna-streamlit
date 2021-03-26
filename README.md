# dna-streamlit

WebApp que faz o alinhamento global par a par de sequências de DNA e analisa a quantidade e a porcentagem de bases e o conteúdo GC das sequências.

Para o alinhamento, foi utilizado o algoritmo de Needleman-Wunsch, considerando apenas match entre as sequências. Cada match foi pontuado +1 na construção da matriz; indels foram considerados como zero.

Para acessar o WebApp, [clique aqui](https://share.streamlit.io/vanleiko/dna-streamlit/main/src/app-dna.py)
