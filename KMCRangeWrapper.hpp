#include <iterator>
#include "./kmc_api/kmc_file.h"

//dummy class to iterate over a kmc database on-disk with iterators

const uint8_t seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

uint64_t pack_kmer(char* kmer, std::size_t k, uint64 mask)
{
	uint64_t pval = 0;
	for (std::size_t i = 0; i < k; ++i) pval = (pval << 2 | seq_nt4_table[(uint8_t)kmer[i]]) & mask;
	return pval;
}

class KmcIterator;
class KMCRangeWrapper
{
    public:
        KMCRangeWrapper(CKMCFile& kmcapi, uint64_t kmask);
        KmcIterator begin() const;//error in BooPHF requires begin() to be const (instead of using cbegin())
        KmcIterator end() const;//error in BooPHF requires end() to be const (instead of using cend())
        ~KMCRangeWrapper();

    protected:
        CKMCFile & kmcdb;
	uint64_t mask;

    private:
        friend class KmcIterator;
};

class KmcIterator
{
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = uint64_t;
        using difference_type   = std::ptrdiff_t;
        using reference         = uint64_t&;
        using pointer           = uint64_t*;
	using iter		= KmcIterator;
        
        explicit KmcIterator(KMCRangeWrapper & kmcrw, std::size_t curr_count, bool is_begin);
        reference operator*();
        pointer operator->();
        iter& operator++();
        iter operator++(int);
        bool operator==(const KmcIterator& rhs) const;
        bool operator!=(const KmcIterator& rhs) const;
        ~KmcIterator();

    private:
        KMCRangeWrapper & parent;
	CKmerAPI kapi;
        std::size_t counter;
        uint64_t current_value;
        uint32_t dummy;
        std::vector<char> buffer;
};

KMCRangeWrapper::KMCRangeWrapper(CKMCFile& kmerapi, uint64_t kmask) : kmcdb(kmerapi), mask(kmask) {}

KmcIterator KMCRangeWrapper::begin() const
{
    kmcdb.RestartListing();
    return KmcIterator(const_cast<KMCRangeWrapper&>(*this), 0, true);
}

KmcIterator KMCRangeWrapper::end() const
{
    return KmcIterator(const_cast<KMCRangeWrapper&>(*this), kmcdb.KmerCount(), false);
}

KMCRangeWrapper::~KMCRangeWrapper() {}

KmcIterator::KmcIterator(KMCRangeWrapper & kmcrw, std::size_t curr_count, bool is_begin) : parent(kmcrw), kapi(kmcrw.kmcdb.KmerLength()), counter(curr_count)
{
    buffer.resize(parent.kmcdb.KmerLength()+1);
    if(is_begin)
    {
        parent.kmcdb.ReadNextKmer(kapi, dummy);
        kapi.to_string(buffer.data());
        current_value = pack_kmer(buffer.data(), parent.kmcdb.KmerLength(), parent.mask);
    }
}

KmcIterator::reference KmcIterator::operator*()
{
//std::cerr << "before dereference access: " << counter << std::endl;
    if(counter >= parent.kmcdb.KmerCount()) throw std::out_of_range("Iterator out of range from dereferencing");
    return current_value;
}

KmcIterator::pointer KmcIterator::operator->()
{
//std::cerr << "before pointer access: " << counter << std::endl;
    if(counter >= parent.kmcdb.KmerCount()) throw std::out_of_range("Iterator out of range from pointer access");
    return &current_value;
}

KmcIterator& KmcIterator::operator++()
{
    parent.kmcdb.ReadNextKmer(kapi, dummy);
    kapi.to_string(buffer.data());
    current_value = pack_kmer(buffer.data(), parent.kmcdb.KmerLength(), parent.mask);
	//std::cerr << "current value = " << current_value << "\n";
	//std::cerr << "before increment: " << counter << " | ";
    ++counter;
	//std::cerr << "after increment: " << counter << std::endl;
    return *this;
}

KmcIterator KmcIterator::operator++(int)
{
    parent.kmcdb.ReadNextKmer(kapi, dummy);
    kapi.to_string(buffer.data());
	//std::cerr << "current value = " << current_value << "\n";    
    current_value = pack_kmer(buffer.data(), parent.kmcdb.KmerLength(), parent.mask);
	//std::cerr << "before increment: " << counter << " | ";
    ++counter;
	//std::cerr << "after increment: " << counter << "\n";
    return *this;
}

bool KmcIterator::operator==(const KmcIterator& rhs) const
{
    return counter == rhs.counter;
}

bool KmcIterator::operator!=(const KmcIterator& rhs) const
{
    return counter != rhs.counter;
}

KmcIterator::~KmcIterator() {}
