#include <iostream>
#include <functional>
#include <vector>
#include <H5Cpp.h>
#include <getopt.h>
#include <assert.h>

static const char OPT_HELP      = 'h';
static const char OPT_REFERENCE = 'r';
static const char OPT_OUTPUT    = 'o';

static const hsize_t CHUNK_SIZE = 32768;

class File : public H5::H5File {
private:
    void init();

protected:
    typedef struct {
	uint8_t  iClass;
	uint64_t iOrderStart;       /* Start index of this class */
	double   fMass;             /* Particle mass */
	double   fSoft;             /* Softening */
	} classEntry;

    H5::CompType m_classType;
    std::vector<classEntry> m_classes;

public:
    File( const char* name, unsigned int flags,
	  const H5::FileCreatPropList& create_plist = H5::FileCreatPropList::DEFAULT,
	  const H5::FileAccPropList& access_plist = H5::FileAccPropList::DEFAULT );
    File( const std::string& name, unsigned int flags,
	  const H5::FileCreatPropList& create_plist = H5::FileCreatPropList::DEFAULT,
	  const H5::FileAccPropList& access_plist = H5::FileAccPropList::DEFAULT );
    File();
    File(const File& original);

    virtual ~File();

    void copyAttributes( const File *src );
    void readAttribute( const std::string &name, double &v );
    void writeAttribute( const std::string &name, double v );


    void readClasses();
    uint64_t getDarkCount() const;
    void fixClasses( const File * RefFile );
    void writeClasses();
    void mergeClasses( const File &splitFile );

    void copyData( const File &src );

    void create( uint64_t N );

public:
    static void doAddAttribute(H5::H5Object&loc, std::string name, void*data);

    // Algorithms
public:

    // compare two File pointers by iOrderStart
class compStart : public std::binary_function<File *,File *,bool> {
    public:
	bool operator()(const File *a, const File *b) const {
	    return a->m_classes[0].iOrderStart < b->m_classes[0].iOrderStart;
	    }
	};

    // Sum particles counts over all files
    class doAddStart {
	uint64_t &m_N;
    public:
	doAddStart( uint64_t &N ) : m_N(N) {
	    m_N = 0;
	    }
	void operator()(const File *a) const {
	    assert( a->m_classes[0].iOrderStart == m_N );
	    m_N += a->getDarkCount();
	    }
	};

    // Read classes from all files
    class doReadClasses {
    public:
	void operator()(File *a) const {
	    a->readClasses();
	    }
	};

    // Fix the classes table from a reference file
    class doFixClasses {
	const File *m_RefFile;
    public:
	doFixClasses( const File * RefFile ) : m_RefFile(RefFile) { }
	void operator()(File *a) const {
	    a->fixClasses( m_RefFile );
	    }
	};

    // Fix the classes table from a reference file
    class doWriteClasses {
    public:
	void operator()(File *a) const {
	    a->writeClasses();
	    }
	};

    // Merge the various class tables
    class doMergeClasses {
	File &m_RefFile;
    public:
	doMergeClasses( File & RefFile ) : m_RefFile(RefFile) { }
	void operator()(File *a) const {
	    m_RefFile.mergeClasses(*a);
	    }
	};

    // Copy the data from one file to another
    class doCopyData {
	File &m_DstFile;
    public:
	doCopyData( File & DstFile ) : m_DstFile(DstFile) { }
	void operator()(File *a) const {
	    m_DstFile.copyData(*a);
	    }
	};

    };
typedef std::vector<File *> FileList;


void File::doAddAttribute(H5::H5Object&loc, std::string name, void*data) {
    double v;
    //hsize_t size;
    H5::Group *d = (H5::Group *)data;

    H5::Attribute a(loc.openAttribute(name));
    //size = loc.getStorageSize();

    //FIXME: double
    a.read( H5::PredType::NATIVE_DOUBLE, &v );

    hsize_t dataSize = 1;
    H5::DataSpace spaceParameters(1,&dataSize,&dataSize);
    H5::Attribute b(d->createAttribute(name,H5::PredType::NATIVE_DOUBLE,
				       spaceParameters));
    b.write( H5::PredType::NATIVE_DOUBLE, &v );
    }

void File::copyAttributes( const File *src ) {
    H5::Group gdst(openGroup( "parameters" ));
    H5::Group gsrc(src->openGroup( "parameters" ));

    gsrc.iterateAttrs(doAddAttribute, NULL, (void *)&gdst);
    }

File::File( const char* name, unsigned int flags,
	    const H5::FileCreatPropList& create_plist,
	    const H5::FileAccPropList& access_plist )
	: H5File(name,flags,create_plist,access_plist),
	m_classType(sizeof(classEntry)) {
    init();
    }

File::File( const std::string& name, unsigned int flags,
	    const H5::FileCreatPropList& create_plist,
	    const H5::FileAccPropList& access_plist)
	: H5File(name,flags,create_plist,access_plist),
	m_classType(sizeof(classEntry)) {
    init();
    }

File::File() : H5File() {
    }


File::File(const File& original) : H5File(original) {
    }

void File::init() {
    m_classType.insertMember("class",HOFFSET(classEntry,iClass),      H5::PredType::NATIVE_UINT8);
    m_classType.insertMember("mass", HOFFSET(classEntry,fMass),       H5::PredType::NATIVE_DOUBLE);
    m_classType.insertMember("soft", HOFFSET(classEntry,fSoft),       H5::PredType::NATIVE_DOUBLE);
    m_classType.insertMember("start",HOFFSET(classEntry,iOrderStart), H5::PredType::NATIVE_UINT64);
    }

File::~File() {
    }

void File::readClasses() {
    hsize_t dims;

    H5::Group group(openGroup( "dark" ));
    H5::DataSet classes(group.openDataSet("classes"));


    assert( classes.getSpace().getSimpleExtentNdims() == 1 );

    classes.getSpace().getSimpleExtentDims (&dims);
    m_classes.insert(m_classes.end(),dims,classEntry());

    classes.read( &(m_classes[0]), m_classType );
    }

void File::writeClasses() {
    H5::Group group(openGroup( "dark" ));
    H5::DataSet classes(group.openDataSet("classes"));
    assert( classes.getSpace().getSimpleExtentNdims() == 1 );

    hsize_t dataSize = m_classes.size();
    hsize_t Zero = 0;
    H5::DataSpace spaceDisk(1,&dataSize);
    H5::DataSpace spaceMem(1,&dataSize);
    classes.extend(&dataSize);
    spaceDisk.selectHyperslab(H5S_SELECT_SET,&dataSize,&Zero);
    spaceMem.selectHyperslab(H5S_SELECT_SET,&dataSize,&Zero);
    classes.write( &(m_classes[0]), m_classType, spaceMem, spaceDisk );
    //classes.flush(H5F_SCOPE_LOCAL);
    }

void File::mergeClasses( const File & splitFile ) {
    if ( m_classes.empty() )
	m_classes.insert(m_classes.end(),
			 splitFile.m_classes.begin(),
			 splitFile.m_classes.end() );
    else {
	std::vector<classEntry>::const_iterator i;
	i = splitFile.m_classes.begin();
	if ( i->iClass == m_classes.back().iClass ) i++;
	m_classes.insert(m_classes.end(), i,
			 splitFile.m_classes.end() );
	}
    }

uint64_t File::getDarkCount() const {
    hsize_t dims[2];

    H5::Group group(openGroup( "dark" ));
    H5::DataSet pos(group.openDataSet("position"));

    assert( pos.getSpace().getSimpleExtentNdims() == 2 );

    pos.getSpace().getSimpleExtentDims (dims);
    assert( dims[1] == 3 );
    return dims[0];
    }

void File::fixClasses( const File *RefFile ) {
    std::vector<classEntry>::const_iterator i,prev;
    uint64_t B, E, X;

    assert( RefFile != NULL );
    assert( !RefFile->m_classes.empty() );

    B = m_classes[0].iOrderStart;
    E = B + getDarkCount();
    for ( i=RefFile->m_classes.begin(); i!=RefFile->m_classes.end(); prev = i, i++ )
	if ( i->iOrderStart > B ) break;
    X = i->iOrderStart;

    m_classes.clear();
    for ( i = prev; i != RefFile->m_classes.end(); i++ ) {
	if ( i->iOrderStart > E ) break;
	m_classes.push_back(classEntry());
	m_classes.back() = *i;
	if ( i->iOrderStart < B )
	    m_classes.back().iOrderStart = B;
	}
    }

void File::readAttribute( const std::string &name, double &v ) {
    H5::Group group(openGroup( "parameters" ));
    H5::Attribute a(group.openAttribute(name));
    a.read( H5::PredType::NATIVE_DOUBLE, &v );
    }

void File::writeAttribute( const std::string &name, double v ) {
    H5::Group group(openGroup( "parameters" ));
    hsize_t dataSize = 1;
    H5::DataSpace spaceParameters(1,&dataSize,&dataSize);
    H5::Attribute a(group.createAttribute(name,H5::PredType::NATIVE_DOUBLE,
					  spaceParameters));
    a.write( H5::PredType::NATIVE_DOUBLE, &v );
    }


void File::create( uint64_t N ) {
    hsize_t chunkSize[2],dataSize[2], maxSize[2];

    maxSize[0] = H5S_UNLIMITED;
    chunkSize[0] = 16;
    dataSize[0] = m_classes.size();
    H5::DataSpace spaceDisk(1,dataSize,maxSize);

    H5::DSetCreatPropList propDisk;
    propDisk.setChunk(1,chunkSize);
    propDisk.setFletcher32();

    H5::Group darkGroup(createGroup( "dark" ));
    H5::DataSet classesDataSet(
	darkGroup.createDataSet("classes", m_classType,
				spaceDisk,propDisk) );

    H5::Group parmGroup(createGroup( "parameters" ));

    maxSize[0] = H5S_UNLIMITED;
    maxSize[1] = 3;
    chunkSize[0] = CHUNK_SIZE;
    chunkSize[1] = 1;
    dataSize[0] = N;
    dataSize[1] = 3;

    H5::DataSpace spaceDisk2(2,dataSize,maxSize);

    H5::DSetCreatPropList propDisk2;
    //propDisk.setLayout(H5D_CHUNKED);
    propDisk2.setChunk(2,chunkSize);
    propDisk2.setFletcher32();

    H5::DataSet positionDataSet(
	darkGroup.createDataSet("position",
				H5::PredType::NATIVE_DOUBLE,
				spaceDisk2,propDisk2) );

    H5::DataSet velocityDataSet(
	darkGroup.createDataSet("velocity",
				H5::PredType::NATIVE_DOUBLE,
				spaceDisk2,propDisk2) );

    }

void File::copyData( const File &src ) {
    uint64_t N = src.getDarkCount();
    uint64_t O = src.m_classes.front().iOrderStart;
    uint64_t i;
    hsize_t n;
    hsize_t off[2], Dims[2];

    H5::Group dstGroup(openGroup( "dark" ));
    H5::Group srcGroup(src.openGroup( "dark" ));

    H5::DataSet dst_pos( dstGroup.openDataSet("position") );
    H5::DataSet dst_vel( dstGroup.openDataSet("velocity") );

    H5::DataSet src_pos( srcGroup.openDataSet("position") );
    H5::DataSet src_vel( srcGroup.openDataSet("velocity") );

    double *buffer = new double[3*CHUNK_SIZE];

    src_pos.getSpace().getSimpleExtentDims(Dims);
    H5::DataSpace src_spaceDisk(2,Dims);
    dst_pos.getSpace().getSimpleExtentDims(Dims);
    H5::DataSpace dst_spaceDisk(2,Dims);


    for ( i=0; i<N; i+=n ) {
	n = N-i > CHUNK_SIZE ? CHUNK_SIZE : N-i;

	Dims[0] = n; Dims[1] = 3;
	H5::DataSpace spaceMem(2,Dims);

	Dims[0] = n; Dims[1] = 3;

	off[0] = off[1] = 0;
	spaceMem.selectHyperslab(H5S_SELECT_SET,Dims,off);

	off[0] = i; off[1] = 0;
	src_spaceDisk.selectHyperslab(H5S_SELECT_SET,Dims,off);

	src_pos.read( buffer, H5::PredType::NATIVE_DOUBLE,
		      spaceMem, src_spaceDisk );

	off[0] = i+O; off[1] = 0;
	dst_spaceDisk.selectHyperslab(H5S_SELECT_SET,Dims,off);

	dst_pos.write( buffer, H5::PredType::NATIVE_DOUBLE,
		       spaceMem, dst_spaceDisk );


	off[0] = i; off[1] = 0;
	src_spaceDisk.selectHyperslab(H5S_SELECT_SET,Dims,off);

	src_vel.read( buffer, H5::PredType::NATIVE_DOUBLE,
		      spaceMem, src_spaceDisk );

	off[0] = i+O; off[1] = 0;
	dst_spaceDisk.selectHyperslab(H5S_SELECT_SET,Dims,off);

	dst_vel.write( buffer, H5::PredType::NATIVE_DOUBLE,
		       spaceMem, dst_spaceDisk );
	}

    delete [] buffer;
    }

int main( int argc, char *argv[] ) {
    bool bHelp = false;
    bool bError = false;
    std::string RefName, OutName;
    File *RefFile = 0;
    FileList files;
    uint64_t n;
    //double dValue;

    //! Parse command line arguments (flags).
    for (;;) {
	int c, option_index=0;

	static struct option long_options[] = {
		{ "help",        0, 0, OPT_HELP
		}, { "reference",   1, 0, OPT_REFERENCE }, { "output",      1, 0, OPT_OUTPUT }, { 0,             0, 0, 0 }
	    };
	c = getopt_long( argc, argv, "hr:o:",
			 long_options, &option_index );
	if ( c == -1 ) break;

	switch (c) {
	case OPT_HELP:
	    bHelp = true;
	    break;
	case OPT_REFERENCE:
	    assert(optarg !=NULL);
	    RefName = optarg;
	    break;
	case OPT_OUTPUT:
	    assert(optarg !=NULL);
	    OutName = optarg;
	    break;
	default:
	    bError = true;
	    }
	}

    if ( bError ) return 1;

    if ( !RefName.empty() ) {
	std::cout << "Opening " << RefName << std::endl;
	RefFile = new File( RefName, H5F_ACC_RDONLY );
	RefFile->readClasses();
	}

    // Open all MDL I/O files
    files.reserve( argc-optind );
    while ( optind < argc ) {
	std::cout << "Opening " << argv[optind] << std::endl;
	files.push_back(
	    new File(argv[optind++],
		     RefName.empty() ? H5F_ACC_RDONLY : H5F_ACC_RDWR));
	}

    // Put the files in iOrder order.
    std::for_each(files.begin(), files.end(), File::doReadClasses());
    std::sort(files.begin(),files.end(), File::compStart());

    // Check that there are no gaps and determine the totals
    std::for_each(files.begin(), files.end(), File::doAddStart(n));
    std::cout << "We have a total of " << n << " dark particles" << std::endl;

    // Now fix the tables if we have a reference file
    if ( RefFile ) {
	std::for_each(files.begin(), files.end(), File::doFixClasses(RefFile));
	std::for_each(files.begin(), files.end(), File::doWriteClasses());
	}

    if ( !OutName.empty() ) {
	std::cout << "Creating " << OutName << std::endl;
	File OutFile( OutName, H5F_ACC_TRUNC );
	OutFile.create(n);
	OutFile.copyAttributes(files.front());
	std::for_each(files.begin(), files.end(), File::doMergeClasses(OutFile));
	OutFile.writeClasses();
	std::for_each(files.begin(), files.end(), File::doCopyData(OutFile));
	}
    }
