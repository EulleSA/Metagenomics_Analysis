
--> É recomendado realizar o processo sem montagem para evitar a formação de quimeras e também para aumentar a cobertura.

--> Interessante falar de shotguns.

--> O KEGG(Kyoto Encyclopedia of Genes and Genomes)[?], 

-->Entretanto, um dos grandes desafios da metagenômica ainda é saber como esses microorganismo interagem entre si.

FERRAMENTAS EXISTENTES

<p align="justify">"Com o passar dos anos diversas ferramentas vêm sendo desenvolvidas com o intuito de auxiliar cientistas e pesquisadores nos seus estudos sobre microbiotas[?]. Exemplos conhecidos como MG-RAST[?] e IMG/G[?], que permitem o usuário fazer as suas análises utilizando o navegador de internet, e o MEGAN6[?]--são muito conhecidos como ferramentas que oferecem a possíbilidade de acesso utilizando a internet,--"</p> 

"Uma forma de entender como esses microorganismos interagem em um determinado meio é análisar as rotas metabólicas que estão presentes." 

"Na metagenômica as duas principais abordagem utilizadas para analisar metagenomas são os baseados em amplicon, onde são buscados genes marcadores de taxonomia no metagenoma, onde os mais conhecidos são o 16s e 18s, e essas sequências são amplificadas para que seja possível a anotação taxonômica[?]."

Possiveis revistas --> BMC Bioinformatics, Oxford Academic(Bioinformatics) e PLOS One(PARECE SER PAGA)

 "A abordagem baseada em shotgun Diferente da abordagem anterior, a baseada em shotgun, onde todo o metagenoma é dividio aleatoriamente em pequenas partes, permitindo assim que seja feita não apenas a anotação taxonômica, mas também a anotação funcional[?]."


"realizando primeiro o alinhamento contra o banco de dados de referência do eggNOG(), e depois a anotação(EXPLICAR MAIS-LER ARTIGO)"

"<p align="justify">O Cactus é construído utilizando o Snakemake, um mecanismo de gerenciamento de workflow que permite a construção de pipeline através de regras definidas em um arquivo chamado 'Snakefile'(Koster and Rahmann 2012). A <span style="color:blue">fig. 1</span> mostra o DIAGRAMA/FLUXOGRAMA da aplicação. Várias ferramentas foram utilizadas neste pipeline: Fastp v0.20.1(Chen et al., 2018) para o controle de qualidade; O eggNOG-mapper v2(Huerta-Cepas et al. ,2017) executa a anotação funcional; Com o output da anotação funcional feita pelo eggNOG-mapper(), é verificado quais reads foram mapeadas para uma ou mais rotas do banco de dados do KEGG orthology(Kanehisa et al., 2017), e então é feita a filtragem das reads do metagenoma pré-processado, mantendo apenas as que possuem rotas metabólicas mapeadas, como mostrado na <span style="color:blue">tabela 1</span>. Para isso, foi utilizada a ferramenta Seqtk(https://github.com/lh3/seqtk); Depois da anotação funcional e da filtragem das reads, partimos para a anotação taxonômica, onde nós buscaremos os organismos pertencentes a esse metagenoma e as rotas metabólicas que eles participam. Primeiro, faz-se o alinhamento  contra o banco de dados NCBI-nr() utilizando a ferramenta DIAMOND v0.9.19; O output gerado pelo DIAMOND é usado no software BASTA(), que utiliza o algoritmo Last Common Ancestral(COLOCAR ALGUM ARTIGO QUE FALE SOBRE ISSO) para fazer a anotação taxonômica.  (FALAR SOBRE A PARTE TAXONOMICA E CRIAR O FLUXOGRAMA/DIAGRAMA DA NOVA METODOLOGIA) </p>"


# <h1 align="justify"> Cactus: A pipeline for identification and visualization of microbial interactions in metabolic pathways using shotgun metagenome sequences. </h1>


## Keywords: metagenomics,pathways,KEGG,shotgun

# <span style="color:yellow"> 1. Introdução </span>

FALAR DE METAGENOMA SHOTGUN

<p align="justify">A metagenômica é uma técnica utilizada para recuperar o material genético de microorganismos diretamente do seu habitat natural, resolvendo um problema anterior no qual os microorganismos precisavam ser cultivados em laboratório para serem estudados[?]. Com o avanço das tecnologias de sequeciamento da nova geração(NGS), principalmente as que realizam sequenciamento de metagenomas baseados em shotgun,onde todo o metagenoma é dividio aleatoriamente em pequenas partes, permitiu-se, através da análise funcional e taxonômica, saber quem são os microorganismos e o que eles fazem em um determinado ambiente(Quince et al., 2017).</p>

ESTUDOS COM ROTAS METABÓLICAS.

<p align="justify">
Um dos grandes desafios da metagenômica é saber como os microorganismos interagem entre si em um dado ambiente. Para compreender isso, diversas ferramentas vêm sendo desenvolvidas para tentar ajudar cientistas e pesquisadores a entender melhor como isso funciona. Entretanto, a maioria delas oferecem a análise funcional e taxonômica separadas, podemos citar como exemplo o MG-RAST(Keegan et al., 2016) e o IMG/G(Chen IA et al., 2016),(que são serviços web e o MEGAN6(Huson DH et al, 2013) que pode ser instalado no seu computador local ou remoto). FALAR ALGO SOBRE AS ROTAS METABÓLICAS -- que oferecem a praticidade e a infraestrutura computacional para que os usuários analisem seus metagenomas, requisitando apenas que o usuário tenha uma boa conexão com a internet e um navegador web para realizar as análises. E também, diferente dessas duas aplicações web, temos o MEGAN6[?] que permite que o usuário instale a ferramenta em seu computador e realize as análises, oferecendo também uma interface de fácil uso.  -- </p> 

SOBRE O MEU PROJETO

Visto isso, nós apresentamos neste artigo Cactus, um pipeline baseado em Snakemake, de fácil uso tanto para pessoas que não tem familiaridade com programação quanto para usuários mais experientes,permitindo a anotação funcional e taxonômica como também a visualização de interações dos organismos presentes nas rotas metabólicas mapeadas na análise de uma determinada microbiota, (a partir do mapeamento das rotas metabólicas mapeadas.) Para isso, foi usado o banco de dados de órtologos do KEGG[], assim como as suas demonstrações de rotas metabólicas. 

# 2. Implemetação

## 2.1 Overview

<p align="justify">O Cactus é construído utilizando Snakemake, um mecanismo de gerenciamento de workflow que permite a construção de pipelines através de regras definidas em um arquivo chamado 'Snakefile'(Koster and Rahmann 2012). O pipeline é dividido em 6 partes: Controle de qualidade, anotação funcional, filtragem das reads e tratamento do output da anotação funcional, alinhamento, anotação taxonômica e geração das redes de interação. Para os passos citados anteriormente foram utilizadas as seguintes ferramentas respectivamente: Fastp v0.20.1(Chen et al., 2018) para o controle de qualidade; eggNOG-mapper v2(Huerta-Cepas et al. ,2017) para a anotação funcional; O tratamento o output da anotação funcional foi feita utilizando a linguagem R() e para a filtragem das reads a ferramenta Seqtk(https://github.com/lh3/seqtk); O DIAMOND v0.9.19 foi usado para realizar o alinhamento contra o banco de dados de referência, NCBI-nr(). O BASTA v?() foi utilizado para realizar a anotação taxonômica; Por fim, a geração das redes de interação foram criadas utilizando a linguagem Python v3 e linguagem R. A <span style="color:blue">fig. 1</span> ilustra os passos e as ferarmentas utilizadas. Mais detalhes serão descritos nos tópicos seguintes.

## 2.2 Controle de qualidade

<p align="justify">Para o controle de qualidade As reads brutas em formato FASTQ foram processadas para realizar a remoção das sequências de baixa qualidade,reads duplicadas e retirada de adaptadores. Para esse processo foi utilizada a ferramenta fastp(CHEN et al., 2018) com PHREAD score menor ou igual a 20.</p>

## 2.3 Anotação Funcional

<p align="justify">A anotação funcional é realizada utilizando a aplicação eggnog-mapper[?], onde recebe como entrada o metagenoma pré-processado, é realizado um alinhamento contra o banco de dados de referência presente no eggnog-mapper, e depois é feito a anotação funcional, dando como arquivo de saída um arquivo tabular no formato ".TSV" e um arquivo de orthologs.</p>

## 2.4 Tratamento dos dados

Utilizando o seqtk(link github), um programa que permite o manuseio e tratamento de arquivos FASTA/FASTQ, é feita a filtragem do metagenomas apenas as reads que foram anotadas para alguma rede metabólica presente no KEGG[?].

## 2.5 Alinhamento 

Com as reads anotadas, utilizando o DIAMOND[?], uma ferramenta tal, que permite alinhar 

## 2.6 Anotação Taxonômica

A anotação taxonômica é feita utilizando o programa BASTA[??], que usa o algoritmo Last Common Ancestral[?] para inferir quais organismos estão presentes no metagenoma.

## 2.7 Geração das redes de interação
Com a anotação funcional e taxonômica realizada, nós conseguimos construir as redes de interação,utilizando KEGGgraph[?], uma biblioteca da linguagem R, presente no repositório do Bioconductor[?], permite a montagem de rotas metabólicas presentes no banco de dados do KEGG em forma de grafos. se lembra de falar de: Utilizando o KEGGgraph, KEGG,script R e Python.

# Metagenoma simulado(dados utilizados)

Para o caso de uso demonstrado a seguir foi gerado um metagenoma simulado composto por xxx sequencias de 7 organismos distintos. Isso só foi possível graças a ferramenta InSilicoSeq, que permite simular dados metagenômicos Illumina. 


# 3. Resultados 






# 4. Discussão


# 5. Conclusão

TENTAR PENSAR NAS DISCUSSÕES E RESULTADOS