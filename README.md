<h1 style="text-align: center;"> Curso de Bioinformática, segundo semestre 2024 </h1>
<h2 style="text-align: center;"> Práctico de RNA-seq </h2>
<p style="text-align: center;"> Jorge Cortés-Miranda & Caren Vega-Retter</p>

---

En este taller exploraremos algunos aspectos fundamentales del análisis de datos provenientes de RNA-seq, iremos desde archivos de secuenciación fastq hasta un análisis de expresión diferencial. Este taller seguirá las recomendaciones expuestas en el artículo de Conesa y colaboradores de 2016 ["A survey for best practices for RNA-seq data analysis"](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8), en partícular se seguirá el caso de un organismo no modelo.

En este taller trabajaremos con datos de secuenciación de *Basilichthys microlepidotus* pejerrey endémico de Chile. Trabajaremos con dos individuos en cuanto a los pasos de mapeo y cuantificación, para luego utilizar una matriz de cuantificación previamente elaborada para el análisis de expresión diferencial.

## **Programas**

Para llevar a cabo el taller utilizaremos programas que debiesen estar instalados en sus equipos. En primer lugar llevaremos a cabo gran parte del taller utilizando **Linux** o **Windows Subsystem for Linux** (WSL). En general basta con que Ubuntu esté instalado como sistema operativo y las herramientas a utilizar dentro del mismo corresponden a **fastqc** para analizar la calidad de las lecturas, **Bowtie2** para llevar a cabo el mapeo de las lecturas a un ensamble *de novo* que se le proporcionará, **samtools** para investigar y filtar los archivos de mapeo y **RSEM** para llevar a cabo el proceso de cuantificación. Por otro lado, es necesario que esté instalado **R studio** y el paquete **DESeq2** para llevar a cabo los análisis de expresión diferencial.

## **Archivos**

Se le proporcionarán 7 archivos, cuatro archivos fastq pertenecientes a dos individuos, uno de San Francisco de Mostazal (SFM1H) y uno de Isla de Maipo (IM2H), y se trabajará con lecturas pareadas o pair-end reads (Por eso son cuatro archivos), un archivo fasta que corresponderá a un ensamble *de novo* de hígado de *Basilichthys microlepidotus*. Se le entregará un archivo que corresponde a una matriz de cuantificación para 6 individuos, 3 que habitan un sitio contaminado en la cuenca del Rio Maipo (MEL) y 3 que habitan un sitio de referencia (SFM). Por último, se le proporcionará un script de R para realizar el análisis de expresión diferencial. 

Cada archivo se encuentra en una carpeta diferente, los archivos .fq se encuentran en samples, mientras que la matriz de cuantificación se encuentra en quantification, el ensamble *de novo* lo encontrarán en transcriptome mientras que en dif_exp encontrarán el script de R para llevar a cabo el análisis de expresión diferencial.

### **Inspección de datos de secuenciación**

Uno de los pasos esenciales cuando se obtienen datos de secuenciación es poder analizar rápidamente la calidad de los mismos, para esto una herramienta como [**fastqc**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) es ideal, para este análisis, en este caso la utilizaremos desde la terminal pero es posible utilizarla en una versión de interfaz gráfica. Para analizar nuestros archivos fastq, iremos a la carpeta samples y ejecutaremos lo siguiente:

```sh
fastqc IM2H_1.fq IM2H_2.fq
```

Esto generará 4 archivos, dos archivos html para cada uno de los archivos fastq analizados y 2 archivos comprimidos en zip. Los archivos html tendrán el reporte de la calidad de las lecturas, mientras que los archivos comprimidos contienen la misma información, pero cada gráfico y datos de forma individual y no integrada como lo sería el caso del archivo html. Si tuvo problemas con este paso puede revisar los archivos que se encuentran en la carpeta quality_res_test. Revisemos los archivos html y veamos los reportes. Discutamos.

Puede repetir el proceso para los archivos provenientes del individuo de SFM.

---

### **Mapeo**

Para llevar a cabo el mapeo utilizaremos [**Bowtie2**](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), el cual es un alineador de secuencias cortas muy útil para datos de RNA-seq. Para utilizarlo y poder mapear nuestras lecturas a una referencia, en este caso nuestro ensamble *de novo* de hígado, es necesario llevar a cabo la creación de una base de datos a partir del archivo fasta del ensamble *de novo*. Para esto nos posicionaremos en nuestra carpeta principal y ejecutaremos lo siguiente en la terminal:

```sh
rsem-prepare-reference --bowtie2 ./transcriptome/hig_transcripts_cdhit80.fasta ./quantification/reference/hig_reference
```
Este programa toma como argumento nuestro ensamble *de novo* en fasta y genera un output con el nombre hig_reference. En realidad, este programa genera la referencia para bowtie2 como también la que utilizaremos para RSEM.

Esto generará una serie de nuevos archivos, los cuales funcionan como índices para cuando se lleva a cabo el mapeo, todos aquellos que terminan en bt2 son utilizados por bowtie2 mientras que el resto servirán a RSEM para la cuantificación.

Si tuvo problemas con este paso se encuentra una carpeta llamada reference_test con todos los archivos que generaría habitualmente en este paso.

Una vez completado este paso, procederemos a realizar el mapeo, para lo cual posicionados en la carpeta RNA_seq ejecutaremos lo siguiente:

```sh
bowtie2 --no-discordant -x ./quantification/reference/hig_reference -1 ./samples/IM2H_1.fq -2 ./samples/IM2H_2.fq -S ./mapping/IM2H.sam 2> ./mapping/mapping_stats.txt
```

Esta línea de comandos nos indica lo siguiente, `--no-discordant` es una opción que evita que se generen alineamiento discordantes, aunque son informados igualmente, `-x` indica el índice de referencia para hacer el mapeo, `-1` y `-2` corresponden a los archivos fastq de las lecturas pareadas para un individuo, `-S` nos indica que el archivo sea de tipo SAM y finalmente `2>` indica que redireccionaremos la salida del stderr hacia un archivo para luego poder inspeccionar las estadísticas del mapeo.

Esta línea de comandos para ejecutar bowtie2 es muy básica y es posible agregar muchas más opciones para lograr el mejor mapeo posible, por lo que es muy importante leer la documentación asociada al momento de buscar la mejor estrategia.

Podemos analizar las estadísticas de nuestro mapeo si nos movemos a la carpeta mapping y utilizando cat revisamos el archivo mapping_stats.txt:

```sh
cat mapping_stats.txt
```
Discutamos lo obtenido.

Puede realizar este paso para el individuo de SFM.

Con esto hemos obtenido el archivo SAM, el cual es un archivo en un formato que puede ser leído por humanos y que da cuenta del mapeo de las lecturas de nuestro individuo.

---

### **Inspección y manipulación de archivo SAM**

Utilizando [**samtools**](http://www.htslib.org/doc/samtools.html) inspeccionaremos nuestro archivo SAM ejecutando lo siguiente:

```sh
samtools view IM2H.sam | head -n 1
```

Esto nos mostrará la primera línea de nuestro archivo SAM, analicémosla.

<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
---


### **Cuantificación**

Es posible llevar a cabo la cuantificación de nuestras lecturas como un proceso separado o gracias a [**RSEM**](https://ycl6.gitbook.io/rna-seq-data-analysis/rna-seq_analysis_workflow/quantification_using_rsem1) como un proceso integrado posterior al mapeo. Gracias a que anteriormente ya creamos una referencia para RSEM, podemos utilizarla y concatenar el mapeo y la cuantificación ejecutando lo siguiente:

```sh
rsem-calculate-expression --bowtie2 -p 4 --paired-end --sort-bam-by-coordinate ./samples/IM2H_1.fq ./samples/IM2H_2.fq ./quantification/reference/hig_reference ./quantification/IM2H_quant
```

Este programa genera el alineamiento con bowtie2, luego llama a samtools y automatiza el proceso de generar archivos BAM que serán utilizados para la cuantificación. Si revisan la carpeta quantification podrán encontrar una serie de archivos, desde los archivos BAM, hasta los resultados de la cuantificación. Ahora abriremos IM2H_quant.genes.results para inspeccionar el archivo.

**Repita este mismo proceso para el individuo de SFM**

Cuando se tiene la cuantificación de varios individuos y deseo generar una matriz que contenga estos datos de cunatifiación para todos ellos, será necesario usar el programa `rsem-generate-data-matrix`, para esto nos moveremos a la carpeta quantification y ejecutaremos lo siguiente:

```sh
rsem-generate-data-matrix SFM1H_quant.genes.results IM2H_quant.genes.results > expression_matrix.tsv
```
podemos explorar este archivo usando head o cat para ver su contenido ¿Qué se ha generado?.

### Expresión Diferencial

En un análisis de expresión diferencial se busca encontrar aquellos genes que, al comparar dos condiciones, por ejemplo, peces que habitan sitios contaminados vs peces que habitan sitios no contaminados, presentan un cambio significativo en su expresión. Para esto es necesario contar con las matrices de expresión que generamos en el paso anterior, donde se consideran a todos los individuos, de ambas condiciones, que deseo comparar. Para este paso utilizaremos una matriz, generada de la misma manera que lo hemos hecho en este taller, pero que contiene 6 individuos y muchas más lecturas. De estos individuos 3 provienen de un sitio contaminado (MEL) y 3 provienen de un sitio de referencia (SFM). Utilizaremos Rstudio, el paquete edgeR y DESeq2 los cuales se encuentran en el repositorio de bioconductor, acá información sobre DESeq2: [https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

El taller continuará con el script de R para el análisis de expresión diferencial donde se detallan los pasos a seguir y se usara el archivo titulado RSEM_count_matrix.txt para dichos pasos.