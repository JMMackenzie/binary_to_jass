#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <numeric>
#include "compress_qmx.h"
#include "compress_integer_elias_gamma_simd.h"
#include <immintrin.h>


uint32_t read_u32(FILE* f)
{
    uint32_t x;
    size_t ret = fread(&x, sizeof(uint32_t), 1, f);
    if (feof(f)) {
        return 0;
    }
    if (ret != 1) {
        exit(EXIT_FAILURE);
    }
    return x;
}

void read_u32s(FILE* f, void* ptr, size_t n)
{
    size_t ret = fread(ptr, sizeof(uint32_t), n, f);
    if (ret != n) {
        exit(EXIT_FAILURE);
    }
}

void read_uint32_list(std::vector<uint32_t>& list, FILE* f)
{
    uint32_t list_len = read_u32(f);
    list.resize(list_len);
    if (list_len != 0) {
        read_u32s(f, list.data(), list_len);
    }
}


//  by D. Lemire, taken from FastPForLib
template<class T>
static void delta(T * data, const size_t size) {
        if (size == 0)
            throw std::runtime_error("delta coding impossible with no value!");
        for (size_t i = size - 1; i > 0; --i) {
            data[i] -= data[i - 1];
        }
}

template<class T>
    static void fastDelta(T * pData, const size_t TotalQty) {
         if (TotalQty < 5) {
             delta(pData, TotalQty);// no SIMD
             return;
         }

         const size_t Qty4 = TotalQty / 4;
         __m128i* pCurr = reinterpret_cast<__m128i*>(pData);
         const __m128i* pEnd = pCurr + Qty4;

         __m128i last = _mm_setzero_si128();
         while (pCurr < pEnd) {
             __m128i a0 = _mm_load_si128(pCurr);
             __m128i a1 = _mm_sub_epi32(a0, _mm_srli_si128(last, 12));
             a1 = _mm_sub_epi32(a1, _mm_slli_si128(a0, 4));
             last = a0;
            
             _mm_store_si128(pCurr++ , a1);
         }

         if (Qty4 * 4 < TotalQty) {
             uint32_t lastVal = _mm_cvtsi128_si32(_mm_srli_si128(last, 12));
             for (size_t i = Qty4 * 4; i < TotalQty; ++i) {
                 uint32_t newVal = pData[i];
                 pData[i] -= lastVal;
                 lastVal = newVal;
             }
         }
    }

 

struct bm25_t {

  static constexpr float b = 0.4;
  static constexpr float k1 = 0.9;
  static constexpr float epsilon_score = 1.0E-6;
  float m_num_docs;
  float m_average_len;

  bm25_t (float num_docs, float average) : m_num_docs(num_docs),
                                           m_average_len(average) {}

  // Computes single-term IDF
  // Assumes term appears once in query (f_qt = 1)
  float query_term_weight(const float f_t) {      // doc frequency
    float idf = std::log((m_num_docs - f_t + 0.5f) / (f_t + 0.5f));
    return std::max(epsilon_score, idf) * (1.0f + k1);
  }

  float doc_term_weight(const float f_t,
                        const float doc_len) {

    float norm_doc_len = doc_len / m_average_len;
    return f_t / (f_t + k1 * (1.0f - b + b * norm_doc_len));
  }

};
constexpr float bm25_t::epsilon_score;


struct posting_t {
  uint32_t docid;
  double score;
  uint32_t q_score;
};

struct plist_t {
  
  std::string term;
  uint32_t term_id;
  std::vector<posting_t> postings_list;

  plist_t(std::string t, uint32_t id) : term(t), term_id(id) {}

};

struct compressed_segment_t {
  uint16_t impact;
  uint32_t frequency;
  uint64_t length;
  uint64_t length_with_padding;
  std::vector<uint32_t> compressed_ids;

  compressed_segment_t (uint16_t imp, uint32_t freq, uint64_t len, uint64_t len_pad, std::vector<uint32_t>&& cids) :
                       impact(imp), frequency(freq), length(len), length_with_padding(len_pad), compressed_ids(cids) {}
};



// Produce CIdoclist.bin
void read_and_write_doclist(std::string docid_file,
                            std::string out_file) {
  std::ifstream docids(docid_file);
  std::ofstream out(out_file);
  std::vector<uint64_t> offsets;
  std::string doc;

  out << "SENTINEL" << '\0';
  offsets.push_back(out.tellp());

  while (docids >> doc) {
    offsets.push_back(out.tellp());
    out << doc << '\0';
  }

  // Now write the offsets
  out.write((char *)&offsets[0], offsets.size() * sizeof(uint64_t));

  uint64_t total_docs = offsets.size();
  // Now write the number of docs
  out.write((char *)&total_docs, sizeof(uint64_t));
}


void score_index(std::string ds2i_prefix,
                 std::string lex_file,
                 const double quant_level,
                 char compression){

	  // Tell JASS that we're using QMX-D1: we need to compute our own deltas
	  ANT_compress *compression_scheme;
    if (compression == 'q')
        compression_scheme = new ANT_compress_qmx();
    else if (compression == 'G')
        compression_scheme = new JASS::compress_integer_elias_gamma_simd();
    else
         compression_scheme = nullptr;
    ANT_compress &compressor = *compression_scheme;


    uint32_t sse_alignment = 16;
    std::vector<plist_t> m_postings_lists;

    // 1. Read lexicon file and start building posting metadata
    std::ifstream lex_data(lex_file);
    std::string term;
    uint32_t term_id = 0;
    while (lex_data >> term) {
      m_postings_lists.emplace_back(term, term_id);
      ++term_id;
    }
    std::cerr << "Read " << m_postings_lists.size() << " terms.\n";

    std::string docs_file = ds2i_prefix + ".docs";
    std::string freqs_file = ds2i_prefix + ".freqs";
    std::string sizes_file = ds2i_prefix + ".sizes";

    auto df = fopen(docs_file.c_str(), "rb");
    auto ff = fopen(freqs_file.c_str(), "rb");
    auto sf = fopen(sizes_file.c_str(), "rb");

    std::vector<uint32_t> doc_list;
    std::vector<uint32_t> freq_list;
    std::vector<uint32_t> size_list;

    double global_maximum = 0.0;

    // Get no. docs    
    read_uint32_list(doc_list,df);
    uint32_t no_docs = doc_list[0];

    // Read sizes
    read_uint32_list(size_list,sf);
  
    // validate sizes
    if (size_list.size() != no_docs) {
      std::cerr << "Something went wrong: Sizes not matched with no. docs." 
                << std::endl;
      exit(EXIT_FAILURE);
    }

    std::cerr << "Processing " << no_docs << " documents." << std::endl;

    // Compute mean document length
    float average_len = std::accumulate(size_list.begin(), size_list.end(), 0.0)/no_docs; 

    bm25_t ranker(float(no_docs), average_len);

    term_id = 0;
    while (!feof(df)) {
      read_uint32_list(doc_list,df);
      read_uint32_list(freq_list,ff);
      size_t doc_n = doc_list.size();
      size_t freq_n = freq_list.size();
      if (doc_n != freq_n) {
        std::cerr << "Files are not aligned.\n";
        exit(EXIT_FAILURE);
      }
      // Exhausted our lists
      if (doc_n == 0) {
        // Ensure validity
        if (term_id != m_postings_lists.size()) {
          std::cerr << "Error: Should have read more data. Exit." << std::endl;
          exit(EXIT_FAILURE);
        }
        break;
      }

      // Query term weight (compute once)
      float w_qt = ranker.query_term_weight(doc_n);

      // Resize our PL
      m_postings_lists[term_id].postings_list.resize(doc_list.size());

      // Score the PL and dump it
      for (size_t i = 0; i < doc_list.size(); ++i) {
        uint32_t doc_id = doc_list[i];
        uint32_t freq = freq_list[i];
        double impact = w_qt * ranker.doc_term_weight(float(freq), float(size_list[doc_id]));
        global_maximum = std::max(global_maximum, impact);
        m_postings_lists[term_id].postings_list[i].docid = doc_id;
        m_postings_lists[term_id].postings_list[i].score = impact;
      }
      ++term_id;
    }

    // We have scored all postings, now we must quantize them
    std::cerr << "Maximum score = " << global_maximum << std::endl;
    std::cerr << "Quantization range = [1, " << quant_level << "]" << std::endl;
    for(size_t i = 0; i < m_postings_lists.size(); ++i) {
      for (size_t j = 0; j < m_postings_lists[i].postings_list.size(); ++j) {
        double score = m_postings_lists[i].postings_list[j].score;
        m_postings_lists[i].postings_list[j].q_score = (uint32_t) std::ceil(quant_level * score / global_maximum);
      }
    }


    // prepare output
    std::ofstream vocab_terms("CIvocab_terms.bin");
    std::ofstream vocab("CIvocab.bin", std::ofstream::binary);
    std::ofstream postings("CIpostings.bin", std::ofstream::binary);
    size_t terms_offset = 0;
    size_t pl_offset = 0;





    postings.write((char *)&compression, sizeof(uint8_t));
    pl_offset = postings.tellp();

    // Now we need to get them out in the expected JASS format
    for(size_t i = 0; i < m_postings_lists.size(); ++i) {

      // Currently processing this list
      std::string current_term = m_postings_lists[i].term;

      // Write the offset in the vocab file
      vocab.write((char *)&terms_offset, sizeof(size_t));
      vocab.write((char *)&pl_offset, sizeof(size_t));

      // Write term to the terms file
      vocab_terms << current_term << '\0';
    
      // Update offsets
      terms_offset = vocab_terms.tellp();

      // Deal with postings
     
      uint32_t max_quant = uint32_t(quant_level);
      std::vector<std::vector<uint32_t>> pl_segments(max_quant);
      for(size_t j = 0; j < m_postings_lists[i].postings_list.size(); ++j) {
        uint32_t score = m_postings_lists[i].postings_list[j].q_score - 1;
        uint32_t docid = m_postings_lists[i].postings_list[j].docid;
        pl_segments[score].push_back(docid + 1); // DocID in JASS is 1->N, 0 is sentinel. Create this offset now.
        // Can now remove this PL...
      }
    
      std::vector<compressed_segment_t> jass_pl;

      // Count non-zero segments and compress them
      uint64_t valid_segments = 0;
      for (size_t j = 0; j < pl_segments.size(); ++j) {
        if (pl_segments[j].size() == 0)
          continue;
      
        // Give the buffer double + some incase of overflow
        std::vector<uint32_t> compressed_segment (2 * pl_segments[j].size() + 1024);
        uint64_t used = compressed_segment.size() * sizeof(compressed_segment[0]);
       
        // Delta encode the document identifiers
        fastDelta(pl_segments[j].data(), pl_segments[j].size()); 
        compressor.encodeArray(pl_segments[j].data(), pl_segments[j].size(), compressed_segment.data(), &used);

        uint64_t used_padded = ((used + sse_alignment - 1) / sse_alignment) * sse_alignment;  // SSE alignment
        // Resize to padded size
        size_t elements = used_padded/sizeof(uint32_t);
        compressed_segment.resize(elements);
        uint16_t impact = j + 1;
        uint32_t frequency = pl_segments[j].size();
        
        jass_pl.emplace_back(impact, frequency, used, used_padded, std::move(compressed_segment));
           
      } 
      // Tell vocab how many segments are here
      valid_segments = jass_pl.size();
      vocab.write((char *)&valid_segments, sizeof(uint64_t));
 
      // We can now write our PL's to file
      // 1. Write the k pointers to headers via arithmetic
      uint32_t header_size = sizeof(uint16_t) + sizeof(uint64_t) + sizeof(uint64_t) + sizeof(uint32_t);
      // Address of the first header is the (number of headers * size of the pointer) + current postings offset.
      uint64_t header_pointer = valid_segments * sizeof(uint64_t) + pl_offset;
      for (size_t i = 0; i < valid_segments; ++i) {
        postings.write((char *)&header_pointer, sizeof(uint64_t));
        header_pointer += header_size;
      }


      uint64_t posting_offset = header_size * (valid_segments + 1) + valid_segments * sizeof(uint64_t) + pl_offset;
      // Ensure SSE alignment is OK
      size_t sse_offset = ((posting_offset + sse_alignment - 1)/sse_alignment) * sse_alignment;
      size_t sse_padding = sse_offset - posting_offset;
      posting_offset = sse_offset;

      // 2. Let's write the header quads now
      for (size_t i = 0; i < jass_pl.size(); ++i) {
        uint16_t impact = jass_pl[i].impact;
        uint32_t frequency = jass_pl[i].frequency;
        uint64_t length = jass_pl[i].length;
        uint64_t length_with_padding = jass_pl[i].length_with_padding;
        uint64_t end = posting_offset + length;
      
        // Write the data
        postings.write((char *)&impact, sizeof(uint16_t));
        postings.write((char *)&posting_offset, sizeof(uint64_t));
        postings.write((char *)&end, sizeof(uint64_t));
        postings.write((char *)&frequency, sizeof(uint32_t));

        // Update offsets for next round
        posting_offset += length_with_padding;

      }

      // Don't forget the zero block
      uint8_t zero_block[] = {0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0};
      postings.write((char *)&zero_block, header_size);

      // More SSE checks
      if (sse_padding != 0) {
        if (sse_padding > sizeof(zero_block)) {
          std::cerr << "Padding is too large (sorry). Fix the code." << std::endl;
          exit(EXIT_FAILURE);
        }
        // Pad some more zeros
        postings.write((char *)zero_block, sse_padding);
      }

      // 4.  Now we write the actual postings lists
      for (size_t i = 0; i < jass_pl.size(); ++i) {
        postings.write((char *)jass_pl[i].compressed_ids.data(), jass_pl[i].length_with_padding);
      }
  
      pl_offset = postings.tellp();

    } 
}

void usage(std::string program) {
  std::cerr << program << " <ds2i_prefix> <docid file> <lexicon file>" 
            << "<max quantized score [256 = 8 bit quant, 512 = 9 bit quant]><q|G [QMX, EliasGammaSIMD]>\n";
  std::cerr << "e.g. " << program << " wsj wsj.documents wsj.terms 256 G\n";
  exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
    if(argc != 6)
      usage(argv[0]);

    // Get input data
    std::string ds2i_prefix = argv[1];
    std::string docid_file = argv[2];
    std::string lexicon_file = argv[3];
    float quantization = std::atof(argv[4]);

    // Do it
    read_and_write_doclist(docid_file, "CIdoclist.bin");
    if (*argv[5] == 'q')
        score_index(ds2i_prefix, lexicon_file, quantization, 'q');
	 else if (*argv[5] == 'G')
        score_index(ds2i_prefix, lexicon_file, quantization, 'G');
    else
        usage(argv[0]);

return (EXIT_SUCCESS);
}

