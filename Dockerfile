# Author: Eric Borenstein

FROM centos:centos7

RUN yum install -y epel-release

RUN yum -y install \
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
	graphviz


ENV APPS_ROOT /apps
RUN mkdir -p ${APPS_ROOT}

###############################################
#BWA = 'bwa/intel/0.7.17'

ENV BWA_VERSION 0.7.17

ENV BWA_HOME ${APPS_ROOT}/bwa/${BWA_VERSION}
ENV PATH ${BWA_HOME}:${PATH}

RUN git clone --branch v${BWA_VERSION} https://github.com/lh3/bwa.git ${BWA_HOME} \
	&& make -C ${BWA_HOME}

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
	&& make \
	&& make install \
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
	&& make \
	&& make install \
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
	&& mv snpEff ${APPS_ROOT}/snpeff/${SNPEFF_VERSION}
				
###############################################
#TRIMMOMATIC = 'trimmomatic/0.36'

ENV TRIMMOMATIC_VERSION 0.36
ENV TRIMMOMATIC_HOME ${APPS_ROOT}/trimmomatic/${TRIMMOMATIC_VERSION}
ENV TRIMMOMATIC_JAR ${TRIMMOMATIC_HOME}/trimmomatic-${TRIMMOMATIC_VERSION}.jar

RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
	&& unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
	&& mkdir -p ${APPS_ROOT}/trimmomatic \
	&& mv Trimmomatic-${TRIMMOMATIC_VERSION} ${TRIMMOMATIC_HOME} \
	&& rm Trimmomatic-${TRIMMOMATIC_VERSION}.zip
			
###############################################
#PYPAIRIX = 'pypairix/intel/0.2.4'

ENV PAIRIX_VERSION 0.2.4
ENV PAIRIX_HOME ${APPS_ROOT/pairix/${PAIRIX_VERSION}
ENV PATH ${PAIRIX_HOME}/bin:${PAIRIX_HOME}/util:${PAIRIX_HOME}/util/bam2pairs:${PATH}

RUN git clone https://github.com/4dn-dcic/pairix --branch ${PAIRIX_VERSION} ${PAIRIX_HOME} \
	&& make -C ${PAIRIX_HOME}

###############################################
#DEEPTOOLS = 'deeptools/3.3.1'

# only python package

###############################################
#JVARKIT = 'jvarkit/base'

ENV JVARKIT_HOME ${APPS_ROOT}/jvarkit
ENV JVARKIT_DIST ${JVARKIT_HOME}/dist

RUN git clone "https://github.com/lindenb/jvarkit.git"

RUN $JVARKIT_HOME/gradlew --project-dir $JVARKIT_HOME sortsamrefname
ENV SORTSAMREFNAME_JAR ${JVARKIT_DIST}/sortsamrefname.jar

RUN $JVARKIT_HOME/gradlew --project-dir $JVARKIT_HOME biostar154220
ENV BIOSTAR_JAR ${JVARKIT_DIST}/biostar154220.jar

###############################################
#PYSAM = 'pysam/intel/python3.6/0.14.1'

# only python package

###############################################
#PILON = 'pilon/1.23'

ENV PILON_VERSION 1.23
ENV PILON_HOME ${APPS_ROOT}/pilon/${PILON_VERSION}
ENV PILON_JAR ${PILON_HOME}/pilon.jar

RUN mkdir -p ${PILON_HOME} \
 && wget https://github.com/broadinstitute/pilon/releases/download/v${PILON_VERSION}/pilon-${PILON_VERSION}.jar -O ${PILON_HOME}/pilon.jar

###############################################

#BCFTOOLS = 'bcftools/intel/1.9'
#BEDTOOLS = 'bedtools/intel/2.27.1'
#NEATGENREADS = 'neat-genreads/v2'
#FREEBAYES = 'freebayes/intel/1.1.0'
#SEQTK = 'seqtk/intel/1.2-r94'
#BIOPYTHON = 'biopython/intel/python3.6/1.72'
#VARSCAN = 'varscan/2.4.2'
#IVAR = 'ivar/1.2.3'
