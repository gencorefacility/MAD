# Author: Eric Borenstein

FROM centos:centos7

RUN  yum update -y \
  && yum install -y epel-release \
  && yum install -y centos-release-scl \
  && yum install -y \
	git \
	wget \
	java-1.8.0-openjdk \
	java-1.8.0-openjdk-devel \
	R \
	autoconf \
	automake \
	make \
	gcc \
	perl-Data-Dumper \
	zlib-devel \
	bzip2 \
	bzip2-devel \
	xz-devel \
	curl-devel \
	openssl-devel \
	ncurses-devel \
	graphviz \
	gsl-devel \ 
	perl-ExtUtils-Embed \
	cmake \
	python2-pip \
	libxml2-devel \
	NLopt* \
	pandoc \
	devtoolset-8-gcc \
	devtoolset-8-gcc-c++ \
   && yum clean all \
   && rm -rf /var/cache/yum

ENV APPS_ROOT /apps
RUN mkdir -p ${APPS_ROOT}

###############################################
# R packages
#RUN R -e "install.packages(c('ggplot2','plyr','tidyverse','ggpubr','MLmetrics','plotrix','rmarkdown'), repos = 'http://cran.us.r-project.org', Ncpus = 6)"

###############################################
#BWA = 'bwa/intel/0.7.17'

ENV BWA_VERSION 0.7.17

ENV BWA_HOME ${APPS_ROOT}/bwa/${BWA_VERSION}
ENV PATH ${BWA_HOME}:${PATH}

RUN git clone --depth 1 --branch v${BWA_VERSION} https://github.com/lh3/bwa.git ${BWA_HOME} \
  && make -j -C ${BWA_HOME}

###############################################
#PICARD = 'picard/2.17.11'

ENV PICARD_VERSION 2.17.11

ENV JAVA_HOME /etc/alternatives/jre
ENV PICARD_HOME ${APPS_ROOT}/picard/${PICARD_VERSION}
ENV PICARD_JAR ${PICARD_HOME}/picard-${PICARD_VERSION}.jar

RUN mkdir -p ${PICARD_HOME}
RUN wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar -O ${PICARD_JAR}

###############################################
#GATK = 'gatk/4.1.7.0'

ENV GATK_VERSION 4.1.7.0

ENV GATK_HOME ${APPS_ROOT}/gatk/${GATK_VERSION}

ENV GATK_LOCAL_JAR ${GATK_HOME}/gatk-package-${GATK_VERSION}-local.jar
ENV GATK_SPARK_JAR ${GATK_HOME}/gatk-package-${GATK_VERSION}-spark.jar
ENV GATK_JAR ${GATK_HOME}/gatk-package-${GATK_VERSION}-local.jar
ENV PATH ${GATK_HOME}:${PATH}

RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip \
  && mkdir ${APPS_ROOT}/gatk \
  && unzip gatk-${GATK_VERSION}.zip \
  && mv gatk-${GATK_VERSION} ${APPS_ROOT}/gatk/${GATK_VERSION} \
  && rm gatk-${GATK_VERSION}.zip

###############################################
#R = 'r/intel/3.4.2'
# INSTALLED MOST CURRENT R

###############################################
#HTSLIB 1.4.1
ENV HTSLIB_VERSION 1.4.1
ENV HTSLIB_HOME ${APPS_ROOT}/htslib/${HTSLIB_VERSION}

ENV MANPATH $MANPATH:${HTSLIB_HOME}/share/man
ENV PATH ${PATH}:${HTSLIB_HOME}/bin
ENV LD_LIBRARY_PATH ${HTSLIB_HOME}/lib:${LD_LIBRARY_PATH}
ENV PKG_CONFIG_PATH ${HTSLIB_HOME}/lib/pkgconfig
ENV HTSLIB_HOME ${HTSLIB_HOME}
ENV HTSLIB_INC ${HTSLIB_HOME}/include
ENV HTSLIB_LIB ${HTSLIB_HOME}/lib

RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
  && tar xjf htslib-${HTSLIB_VERSION}.tar.bz2 \
  && rm htslib-${HTSLIB_VERSION}.tar.bz2 \
  && cd htslib-${HTSLIB_VERSION} \
  && autoheader \
  && autoconf  \
  && ./configure --prefix=${HTSLIB_HOME} \
  && make -j \
  && make -j install \
  && cd / \
  && rm -Rf /htslib-${HTSLIB_VERSION}

###############################################
#SAMTOOLS = 'samtools/intel/1.9'

ENV SAMTOOLS_VERSION 1.9
ENV SAMTOOLS_HOME ${APPS_ROOT}/samtools/${SAMTOOLS_VERSION}

ENV MANPATH ${SAMTOOLS_HOME}/share/man
ENV PATH ${SAMTOOLS_HOME}/bin:${PATH}
ENV LD_LIBRARY_PATH ${SAMTOOLS_HOME}/lib:${LD_LIBRARY_PATH}
ENV SAMTOOLS_HOME ${SAMTOOLS_HOME}
ENV SAMTOOLS_INC ${SAMTOOLS_HOME}/include
ENV SAMTOOLS_LIB ${SAMTOOLS_HOME}/lib

RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && cd samtools-${SAMTOOLS_VERSION} \
  && autoheader \
  && autoconf -Wno-syntax \
  && ./configure --prefix=${SAMTOOLS_HOME} \
  && make -j \
  && make -j install \
  && cd / \
  && rm -Rf samtools-${SAMTOOLS_VERSION}
	
###############################################
#SNPEFF = 'snpeff/4.3i'

ENV SNPEFF_VERSION 4_3i
ENV SNPEFF_HOME ${APPS_ROOT}/snpeff/${SNPEFF_VERSION}

ENV SNPEFF_JAR ${SNPEFF_HOME}/snpEff.jar
ENV SNPSIFT_JAR ${SNPEFF_HOME}/SnpSift.jar

RUN wget -O snpEff_v${SNPEFF_VERSION}_core.zip  https://sourceforge.net/projects/snpeff/files/snpEff_v${SNPEFF_VERSION}_core.zip/download# \
  && mkdir ${APPS_ROOT}/snpeff \
  && unzip snpEff_v${SNPEFF_VERSION}_core.zip \
  && mv snpEff ${APPS_ROOT}/snpeff/${SNPEFF_VERSION} \
	&& rm snpEff_v${SNPEFF_VERSION}_core.zip
				
###############################################
#TRIMMOMATIC = 'trimmomatic/0.36'

ENV TRIMMOMATIC_VERSION 0.36
ENV TRIMMOMATIC_HOME ${APPS_ROOT}/trimmomatic/${TRIMMOMATIC_VERSION}
ENV TRIMMOMATIC_JAR ${TRIMMOMATIC_HOME}/trimmomatic-${TRIMMOMATIC_VERSION}.jar

RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
  && unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
  && mkdir -p ${APPS_ROOT}/trimmomatic \
  && mv Trimmomatic-${TRIMMOMATIC_VERSION} ${TRIMMOMATIC_HOME} \
	&& rm Trimmomatic-${TRIMMOMATIC_VERSION}.zip
	
###############################################
#PYPAIRIX = 'pypairix/intel/0.2.4'

ENV PAIRIX_VERSION 0.2.4
ENV PAIRIX_HOME ${APPS_ROOT}/pairix/${PAIRIX_VERSION}
ENV PATH ${PAIRIX_HOME}/bin:${PAIRIX_HOME}/util:${PAIRIX_HOME}/util/bam2pairs:${PATH}

RUN git clone --depth 1 https://github.com/4dn-dcic/pairix --branch ${PAIRIX_VERSION} ${PAIRIX_HOME} \
  && make -j -C ${PAIRIX_HOME}

###############################################
#JVARKIT = 'jvarkit/base'

ENV JVARKIT_HOME ${APPS_ROOT}/jvarkit
ENV JVARKIT_DIST ${JVARKIT_HOME}/dist

RUN git clone --depth 1 "https://github.com/lindenb/jvarkit.git" ${JVARKIT_HOME}

RUN $JVARKIT_HOME/gradlew --project-dir $JVARKIT_HOME sortsamrefname
ENV SORTSAMREFNAME_JAR ${JVARKIT_DIST}/sortsamrefname.jar

RUN $JVARKIT_HOME/gradlew --project-dir $JVARKIT_HOME biostar154220
ENV BIOSTAR_JAR ${JVARKIT_DIST}/biostar154220.jar

###############################################
#PILON = 'pilon/1.23'

ENV PILON_VERSION 1.23
ENV PILON_HOME ${APPS_ROOT}/pilon/${PILON_VERSION}
ENV PILON_JAR ${PILON_HOME}/pilon.jar

RUN mkdir -p ${PILON_HOME} \
 && wget https://github.com/broadinstitute/pilon/releases/download/v${PILON_VERSION}/pilon-${PILON_VERSION}.jar -O ${PILON_HOME}/pilon.jar

###############################################
#BCFTOOLS = 'bcftools/intel/1.9'

ENV BCFTOOLS_VERSION 1.9
ENV BCFTOOLS_HOME ${APPS_ROOT}/bcftools/${BCFTOOLS_VERSION}
ENV BCFTOOLS_PLUGINS ${BCFTOOLS_HOME}/plugins
ENV PATH ${BCFTOOLS_HOME}:${PATH}
ENV MANPATH ${BCFTOOLS_HOME}/share/man:${MANPATH}
ENV LD_LIBRARY_PATH ${BCFTOOLS_HOME}/lib:${LD_LIBRARY_PATH}

RUN git clone --depth 1 git://github.com/samtools/bcftools.git --branch ${BCFTOOLS_VERSION} ${BCFTOOLS_HOME} \
  && cd ${BCFTOOLS_HOME} \
  && autoheader && autoconf && ./configure --with-htslib=${HTSLIB_HOME} \
  && make -j

###############################################
#BEDTOOLS = 'bedtools/intel/2.27.1'

ENV BEDTOOLS_VERSION 2.27.1
ENV BEDTOOLS_HOME ${APPS_ROOT}/bedtools/${BEDTOOLS_VERSION}
ENV PATH ${BEDTOOLS_HOME}/bin:${PATH}

RUN wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz \
  && tar -zxvf bedtools-${BEDTOOLS_VERSION}.tar.gz \
  && make -j -C bedtools2 \
  && mkdir -p ${BEDTOOLS_HOME}/bin \
  && mv bedtools2/bin/* ${BEDTOOLS_HOME}/bin \
  && rm -Rf bedtools2 bedtools-${BEDTOOLS_VERSION}.tar.gz

###############################################
#SEQTK = 'seqtk/intel/1.2-r94'

ENV SEQTK_VERSION 1.2
ENV SEQTK_HOME ${APPS_ROOT}/seqtk/${SEQTK_VERSION}
ENV PATH ${SEQTK_HOME}:${PATH}

RUN git clone --depth 1 https://github.com/lh3/seqtk.git --branch v${SEQTK_VERSION} ${SEQTK_HOME} \
  && make -j -C ${SEQTK_HOME}

###############################################
#FREEBAYES = 'freebayes/intel/1.1.0'

ENV FREEBAYES_VERSION 1.1.0
ENV FREEBAYES_HOME ${APPS_ROOT}/freebayes/${FREEBAYES_VERSION}
ENV PATH ${FREEBAYES_HOME}/bin:${PATH}

RUN git clone --depth 1 --branch v${FREEBAYES_VERSION} --recursive git://github.com/ekg/freebayes.git ${FREEBAYES_HOME} \
  && make -C ${FREEBAYES_HOME}
  
###############################################
#VARSCAN = 'varscan/2.4.2'

ENV VARSCAN_VERSION 2.4.2
ENV VARSCAN_HOME ${APPS_ROOT}/varscan/${VARSCAN_VERSION}
ENV VARSCAN_JAR ${VARSCAN_HOME}/VarScan.v${VARSCAN_VERSION}.jar

RUN mkdir -p ${VARSCAN_HOME} \
  && wget https://github.com/dkoboldt/varscan/releases/download/${VARSCAN_VERSION}/VarScan.v${VARSCAN_VERSION}.jar -O ${VARSCAN_JAR}


###############################################
#IVAR = 'ivar/1.2.3'

ENV IVAR_VERSION 1.2.3
ENV IVAR_HOME ${APPS_ROOT}/ivar/${IVAR_VERSION}

RUN git clone --depth 1 https://github.com/andersen-lab/ivar.git --branch v${IVAR_VERSION} ${IVAR_HOME} \
  && cd ${IVAR_HOME} \
  && ./autogen.sh \
  && ./configure --with-hts=${HTSLIB_HOME} \
  && scl enable devtoolset-8 -- make -j \
  && scl enable devtoolset-8 -- make install \
  && cd /

###############################################
#NEATGENREADS = 'neat-genreads/v2'
#Requires: Python 2.7, Numpy 1.9.1+

ENV NEATGENREADS_VERSION 2.0
ENV NEATGENREADS_HOME ${APPS_ROOT}/neat-genreads/${NEATGENREADS_VERSION}

RUN pip2 install --upgrade pip==20.3.1 \
  && pip2 install numpy==1.16.6 \
  && git clone --depth 1 --branch v${NEATGENREADS_VERSION} https://github.com/zstephens/neat-genreads.git ${NEATGENREADS_HOME} \
  && echo 'alias genReads.py="python2 ${NEATGENREADS_HOME}/genReads.py"' >> /etc/profile.d/neatgen.sh \
  && echo 'alias mergeJobs.py="python2 ${NEATGENREADS_HOME}/mergeJobs.py"' >> /etc/profile.d/neatgen.sh \
  && echo 'alias computeFraglen.py="python2 ${NEATGENREADS_HOME}/utilities/computeFraglen.py"' >> /etc/profile.d/neatgen.sh \
  && echo 'alias computeGC.py="python2 ${NEATGENREADS_HOME}/utilities/computeGC.py"' >> /etc/profile.d/neatgen.sh \
  && echo 'alias genMutModel.py="python2 ${NEATGENREADS_HOME}/utilities/genMutModel.py"' >> /etc/profile.d/neatgen.sh \
  && echo 'alias genSeqErrorModel.py="python2 ${NEATGENREADS_HOME}/utilities/genSeqErrorModel.py"' >> /etc/profile.d/neatgen.sh \
  && echo 'alias plotMutModel.py="python2 ${NEATGENREADS_HOME}/utilities/plotMutModel.py"' >> /etc/profile.d/neatgen.sh \
  && echo 'alias validateBam.py="python2 ${NEATGENREADS_HOME}/utilities/validateBam.py"' >> /etc/profile.d/neatgen.sh \
  && echo 'alias validateFQ.py="python2 ${NEATGENREADS_HOME}/utilities/validateFQ.py"' >> /etc/profile.d/neatgen.sh \
  && echo 'alias vcf_compare_OLD.py="python2 ${NEATGENREADS_HOME}/utilities/vcf_compare_OLD.py"' >> /etc/profile.d/neatgen.sh

###############################################
#CONDA + 
#biopython==1.72 deeptools==3.3.1 pysam==0.14.1
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh \
  && bash ~/miniconda.sh -b -p ${APPS_ROOT}/miniconda \
  && eval "$(${APPS_ROOT}/miniconda/bin/conda shell.bash hook)" \
  && conda config --add channels defaults \
  && conda config --add channels bioconda \
  && conda config --add channels conda-forge \
  && conda install biopython==1.72 deeptools==3.3.1 pysam==0.14.1 python==3.6 \
  && echo '. ${APPS_ROOT}/miniconda/etc/profile.d/conda.sh' >> /etc/profile.d/miniconda.sh \
  && conda clean -t -y
		
###############################################
#PYSAM = pysam/intel/0.10.0

RUN pip2 install pysam==0.10.0
RUN echo 'conda activate base' >> /etc/profile.d/miniconda.sh

###############################################
#lofreq_star/2.1.3.1
ENV LOFREQ_STAR_VERSION 2.1.3.1
ENV LOFREQ_STAR_HOME ${APPS_ROOT}/lofreq_star/${LOFREQ_STAR_VERSION}
ENV PATH ${PATH}:${LOFREQ_STAR_HOME}/bin

RUN wget https://github.com/CSB5/lofreq/raw/master/dist/lofreq_star-${LOFREQ_STAR_VERSION}_linux-x86-64.tgz \
  && tar -zxvf lofreq_star-${LOFREQ_STAR_VERSION}_linux-x86-64.tgz \
	&& mkdir -p ${LOFREQ_STAR_HOME} \
	&& mv ./lofreq_star-2.1.3.1/bin  ${LOFREQ_STAR_HOME}/ 
	
ENV PATH ${PATH}:${APPS_ROOT}/miniconda/bin/

RUN R -e "install.packages(c('ggplot2','plyr','tidyverse','ggpubr','MLmetrics','plotrix','rmarkdown'), repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('devtools', repos = 'http://cran.us.r-project.org')"
RUN R -e "require(devtools); install_version('tidyr', version = '1.0.0', repos = 'http://cran.us.r-project.org')"
RUN R -e "require(devtools); install_version('forcats', version = '0.4.0', repos = 'http://cran.us.r-project.org')"
RUN R -e "require(devtools); install_version('readr', version = '1.3.1', repos = 'http://cran.us.r-project.org')"
RUN R -e "require(devtools); install_version('plyr', version = '1.8.5', repos = 'http://cran.us.r-project.org')"
RUN R -e "require(devtools); install_version('purrr', version = '0.3.3', repos = 'http://cran.us.r-project.org')"

RUN mkdir ${APP_ROOT}/scripts
ADD bin/analyze_af_report.Rmd ${APPS_ROOT}/scripts/analyze_af_report.Rmd

#RUN R -e "require(devtools); install_version('ggplot2', version = '3.2.1', repos = 'http://cran.us.r-project.org')"
#RUN R -e "require(devtools); install_version('jsonlite', version = '1.6', repos = 'http://cran.us.r-project.org')"
#RUN R -e "require(devtools); install_version('tidyverse', version = '1.2.1', repos = 'http://cran.us.r-project.org')"

RUN R -e "install.packages(c('tinytex'), repos = 'http://cran.us.r-project.org', Ncpus = 6)"
RUN R -e "tinytex::install_tinytex()"


