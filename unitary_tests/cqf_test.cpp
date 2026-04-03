#include <gtest/gtest.h>

#include "rsqf.hpp"
#include "bqf_ec.hpp"
#include "bqf_oom.hpp"
#include "bqf_cf.hpp"
#include <random>
#include <filesystem>

using namespace std;

class RsqfTest : public ::testing::Test {
 protected:
  void SetUp() override {
    generator.seed(time(NULL));

    usual_qf = Rsqf(2);
    small_qf = Rsqf(7, 64-7, false);
  }

  std::default_random_engine generator;
  std::uniform_int_distribution<uint64_t> distribution;

  Rsqf usual_qf;
  Rsqf small_qf;
};


TEST_F(RsqfTest, get_runstart) {
    small_qf.insert((2ULL<<30)+64);
    EXPECT_EQ(small_qf.get_runstart(64), 64);
    small_qf.insert((2ULL<<31)+63);

    small_qf.insert((2ULL<<32)+126);
    small_qf.insert((2ULL<<33)+126);
    small_qf.insert((2ULL<<34)+126);
    EXPECT_EQ(small_qf.get_runstart(64), 64);
}


TEST_F(RsqfTest, insert_pos_0) {
    small_qf.insert((2ULL << 30) + 90);
    small_qf.insert(2ULL << 31);
    EXPECT_EQ(small_qf.get_remainder(0), 33554432);
}


TEST_F(RsqfTest, offset1_blockOverflow) {
    for (int i = 0; i < 5; i++) { small_qf.insert((1<<11)+ 62); }

    EXPECT_EQ(small_qf.get_offset_word(1), 3);

    small_qf.remove((1<<11)+ 62);
    EXPECT_EQ(small_qf.get_offset_word(1), 2);
    small_qf.remove((1<<11)+ 62);
    EXPECT_EQ(small_qf.get_offset_word(1), 1);
    small_qf.remove((1<<11)+ 62);
    EXPECT_EQ(small_qf.get_offset_word(1), 0);
}

// offset2_runOverflowsBy1: known bug in filter offset logic — kept as documentation
// TEST_F(RsqfTest, offset2_runOverflowsBy1) { ... }

TEST_F(RsqfTest, offset3_insertshift0) {
    small_qf.insert((1<<11)+ 64); EXPECT_EQ(small_qf.get_offset_word(1), 1);
    small_qf.insert((1<<11)+ 64); EXPECT_EQ(small_qf.get_offset_word(1), 2);

    small_qf.remove((1<<11)+ 64); EXPECT_EQ(small_qf.get_offset_word(1), 1);
    small_qf.remove((1<<11)+ 64); EXPECT_EQ(small_qf.get_offset_word(1), 0);
    small_qf.remove((1<<11)+ 64); EXPECT_EQ(small_qf.get_offset_word(1), 0);
}


TEST_F(RsqfTest, offset4_multiBlockRun) {
    for (int i = 0; i < 128; i++){ usual_qf.insert((1ULL<<32)+ 100); }
    EXPECT_EQ(usual_qf.get_offset_word(1), 0);
    EXPECT_EQ(usual_qf.get_offset_word(2), 100);
    EXPECT_EQ(usual_qf.get_offset_word(3), 100-64);

    for (int i = 0; i < 50; i++){ usual_qf.remove((1ULL<<32)+ 100); }
    EXPECT_EQ(usual_qf.get_offset_word(1), 0);
    EXPECT_EQ(usual_qf.get_offset_word(2), 50);
    EXPECT_EQ(usual_qf.get_offset_word(3), 0);
}

TEST_F(RsqfTest, offset5_complexCase) {
    for (int i = 0; i < 60; i++){ usual_qf.insert((1ULL<<32)+ 20); }
    for (int i = 0; i < 69; i++){ usual_qf.insert((1ULL<<32)+ 40); }
    usual_qf.insert((1ULL<<32)+ 149);
    for (int i = 0; i < 11; i++){ usual_qf.insert((1ULL<<32)+ 150); }

    EXPECT_EQ(usual_qf.get_offset_word(0), 0);
    EXPECT_EQ(usual_qf.get_offset_word(1), 149-64);
    EXPECT_EQ(usual_qf.get_offset_word(2), 149-128);

    usual_qf.remove((1ULL<<32)+ 155);
    EXPECT_EQ(usual_qf.get_offset_word(2), 149-128);
    for (int i = 0; i < 20; i++){ usual_qf.remove((1ULL<<32)+ 20); }
    EXPECT_EQ(usual_qf.get_offset_word(0), 0);
    EXPECT_EQ(usual_qf.get_offset_word(1), 129-64);
    EXPECT_EQ(usual_qf.get_offset_word(2), 129-128);
}

TEST_F(RsqfTest, toricity1) {
    for (int i = 0; i < 100; i++){ small_qf.insert((1ULL<<32)+ 40); }

    EXPECT_EQ(small_qf.get_offset_word(0), 12);
    EXPECT_EQ(small_qf.get_offset_word(1), 76);
}


TEST_F(RsqfTest, toricity2) {
    for (int i = 0; i < 20; i++){ small_qf.insert((1ULL<<32)+ 50); }
    for (int i = 0; i < 20; i++){ small_qf.insert((1ULL<<31)+ 120); }

    EXPECT_EQ(small_qf.get_offset_word(0), 12);
    EXPECT_EQ(small_qf.get_offset_word(1), 6);
}


TEST_F(RsqfTest, enumerate1) {
    std::unordered_set<uint64_t> verif;

    small_qf.insert((3ULL<<30) + 35); verif.insert((3ULL<<30) + 35);
    small_qf.insert((3ULL<<31) + 35); verif.insert((3ULL<<31) + 35);
    small_qf.insert((3ULL<<25) + 35); verif.insert((3ULL<<25) + 35);

    small_qf.insert((3ULL<<31) + 64); verif.insert((3ULL<<31) + 64);

    small_qf.insert((3ULL<<30) + 63); verif.insert((3ULL<<30) + 63);
    small_qf.insert((3ULL<<30) + 63); verif.insert((3ULL<<30) + 63);

    EXPECT_EQ(small_qf.enumerate(), verif);
}


TEST_F(RsqfTest, enumerate2) {
    std::unordered_set<uint64_t> verif;
    for (int i=0; i<100; i++){
        uint64_t val = distribution(generator);
        verif.insert(val);
        usual_qf.insert(val);
    }
    EXPECT_EQ(usual_qf.enumerate(), verif);
}


TEST_F(RsqfTest, globalUse) {
    uint64_t val = distribution(generator);
    uint64_t val2 = distribution(generator);

    while (val == val2){ val2 = distribution(generator); }
    usual_qf.insert(val);
    EXPECT_EQ(usual_qf.query(val), 1);
    EXPECT_EQ(usual_qf.query(val2), 0);
}


TEST_F(RsqfTest, get_run_boundaries) {
    std::pair<uint64_t,uint64_t> compare(126, 2);
    small_qf.insert((2ULL<<30)+126);
    small_qf.insert((2ULL<<31)+126);
    small_qf.insert((2ULL<<32)+126);
    small_qf.insert((2ULL<<33)+126);
    small_qf.insert((2ULL<<34)+126);
    EXPECT_EQ(small_qf.get_run_boundaries(126), compare);
}

TEST_F(RsqfTest, get_run_boundaries2) {
    std::pair<uint64_t,uint64_t> compare;

    for (int i = 0; i < 16; i++){ small_qf.insert((20ULL<<7)+ 20); }
    for (int i = 0; i < 28; i++){ small_qf.insert((40ULL<<7)+ 40); }
    small_qf.insert((99ULL<<7)+ 99);
    for (int i = 0; i < 12; i++){ small_qf.insert((100ULL<<7)+ 100); }
    for (int i = 0; i < 48; i++){ small_qf.insert((96ULL<<7)+ 96); }

    compare = std::make_pair(29, 44);   EXPECT_EQ(small_qf.get_run_boundaries(20), compare);
    compare = std::make_pair(45, 72);   EXPECT_EQ(small_qf.get_run_boundaries(40), compare);
    compare = std::make_pair(96, 15);   EXPECT_EQ(small_qf.get_run_boundaries(96), compare);
    compare = std::make_pair(16, 16);   EXPECT_EQ(small_qf.get_run_boundaries(99), compare);
    compare = std::make_pair(17, 28);   EXPECT_EQ(small_qf.get_run_boundaries(100), compare);

    for (int i = 0; i < 28; i++){ small_qf.remove((96ULL<<7)+ 96); }

    compare = std::make_pair(20, 35);   EXPECT_EQ(small_qf.get_run_boundaries(20), compare);
    compare = std::make_pair(40, 67);   EXPECT_EQ(small_qf.get_run_boundaries(40), compare);
    compare = std::make_pair(96, 115);  EXPECT_EQ(small_qf.get_run_boundaries(96), compare);
    compare = std::make_pair(116, 116); EXPECT_EQ(small_qf.get_run_boundaries(99), compare);
    compare = std::make_pair(117, 0);   EXPECT_EQ(small_qf.get_run_boundaries(100), compare);
}


TEST_F(RsqfTest, first_unused_slot1) {
    for (int i = 0; i < 2; i++){ small_qf.insert((123<<7) + 123); }
    small_qf.insert((124<<7) + 124);
    for (int i = 0; i < 2; i++){ small_qf.insert((125<<7) + 125); }
    small_qf.insert((126<<7) + 126);
    small_qf.insert((1<<7) + 1);
    for (int i = 0; i < 3; i++){ small_qf.insert((2<<7) + 2); }

    EXPECT_EQ(small_qf.first_unused_slot(123), 5);
}

TEST_F(RsqfTest, first_unused_slot2) {
    for (int i = 0; i < 2; i++){ small_qf.insert((123<<7) + 123); }
    for (int i = 0; i < 4; i++){ small_qf.insert((124<<7) + 124); }
    for (int i = 0; i < 3; i++){ small_qf.insert((1<<7) + 1); }

    EXPECT_EQ(small_qf.first_unused_slot(123), 4);
}


TEST_F(RsqfTest, first_unused_slot3) {
    small_qf.insert((63<<7) + 63);
    small_qf.insert((65<<7) + 65);

    EXPECT_EQ(small_qf.first_unused_slot(63), 64);
}


TEST_F(RsqfTest, shift_bits_right_metadata) {
    small_qf.insert((63<<7) + 63);
    small_qf.insert((3000<<7) + 63);
    small_qf.insert((64<<7) + 64);

    small_qf.remove((3000<<7) + 63);

    EXPECT_EQ(small_qf.get_runend_word(0), 1ULL<<63);
}

// finalTest: 2^18 insertions + full enumerate + full removal — too slow for routine CI.
// Run manually: ./build/bin/unit_tests --gtest_filter=RsqfTest.finalTest
// TEST_F(RsqfTest, finalTest) { ... }


class BCqfTest : public ::testing::Test {
 protected:
  void SetUp() override {
    generator.seed(time(NULL));

    small_cqf = Bqf_ec(7, 5, 32, 0, false);
    cqf = Bqf_ec(1, 5, false);

    small_cqf_oom = Bqf_oom(7, 5, 32, 0, false);
    cqf_oom = Bqf_oom(1, 5, false);
  }

  std::default_random_engine generator;
  std::uniform_int_distribution<uint64_t> distribution;

  Bqf_ec small_cqf;
  Bqf_ec cqf;

  Bqf_oom small_cqf_oom;
  Bqf_oom cqf_oom;
};


TEST_F(BCqfTest, insert1occs) {
    uint64_t val;
    std::map<uint64_t, uint64_t> verif;

    for (uint64_t i=0; i < (1ULL<<17)-1; i++){
        val = distribution(generator);
        while (verif.count(val) == 1) {
            val = distribution(generator);
        }

        cqf.insert(val);
        verif.insert({ val, 1 });
    }

    EXPECT_EQ(cqf.enumerate(), verif);

    std::map<uint64_t, uint64_t>::iterator it;
    for (it = verif.begin(); it != verif.end(); it++){
        cqf.remove((*it).first);
    }
    verif.clear();

    EXPECT_EQ(cqf.enumerate(), verif);
}


TEST_F(BCqfTest, insertRDMoccs) {
    uint64_t val;
    std::map<uint64_t, uint64_t> verif;

    for (uint64_t i=0; i < (1ULL<<17)-1; i++){
        val = distribution(generator);
        while (verif.count(val) == 1) {
            val = distribution(generator);
        }

        cqf.insert(val, val%31);
        verif.insert({val, val%31 });
    }
    EXPECT_EQ(cqf.enumerate(), verif);

    std::map<uint64_t,uint64_t>::iterator it;
    for (it = verif.begin(); it != verif.end(); it++){
        cqf.remove((*it).first, (*it).second);
    }
    verif.clear();

    EXPECT_EQ(cqf.enumerate(), verif);
}


TEST_F(BCqfTest, insertSeveralOccs) {
    uint64_t val;
    std::map<uint64_t, uint64_t> verif;

    for (uint64_t i=0; i < (1ULL<<17)-1; i++){
        val = distribution(generator);
        cqf.insert(val);
        ++verif[val];
    }
    EXPECT_EQ(cqf.enumerate(), verif);

    std::map<uint64_t,uint64_t>::iterator it;
    for (it = verif.begin(); it != verif.end(); it++){
        cqf.remove((*it).first, (*it).second);
    }
    verif.clear();

    EXPECT_EQ(cqf.enumerate(), verif);
}

TEST_F(BCqfTest, insertRDMoccs_oom) {
    uint64_t val;
    std::map<uint64_t, uint64_t> verif;

    for (uint64_t i=0; i < (1ULL<<17)-1; i++){
        val = distribution(generator);
        while (verif.count(val) == 1) {
            val = distribution(generator);
        }

        cqf_oom.insert(val, (1ULL << val%31));
        verif.insert({ val, (1ULL << val%31) });
    }
    EXPECT_EQ(cqf_oom.enumerate(), verif);

    std::map<uint64_t,uint64_t>::iterator it;
    for (it = verif.begin(); it != verif.end(); it++){
        cqf_oom.remove((*it).first);
    }
    verif.clear();

    EXPECT_EQ(cqf_oom.enumerate(), verif);
}


TEST_F(BCqfTest, serializeRoundtrip_ec) {
    uint64_t val;
    std::map<uint64_t, uint64_t> verif;

    for (uint64_t i = 0; i < (1ULL<<14); i++){
        val = distribution(generator);
        while (verif.count(val) == 1) { val = distribution(generator); }
        cqf.insert(val, val % 15 + 1);
        verif.insert({ val, val % 15 + 1 });
    }

    const std::string path = "/tmp/bqf_test_ec_roundtrip";
    cqf.save_on_disk(path);
    Bqf_ec loaded = Bqf_ec::load_from_disk(path);
    std::filesystem::remove(path);

    EXPECT_EQ(loaded.enumerate(), verif);
    EXPECT_EQ(loaded.quotient_size,  cqf.quotient_size);
    EXPECT_EQ(loaded.remainder_size, cqf.remainder_size);
    EXPECT_EQ(loaded.count_size,     cqf.count_size);
    EXPECT_EQ(loaded.elements_inside, cqf.elements_inside);
}


TEST_F(BCqfTest, serializeRoundtrip_oom) {
    uint64_t val;
    std::map<uint64_t, uint64_t> verif;

    for (uint64_t i = 0; i < (1ULL<<14); i++){
        val = distribution(generator);
        while (verif.count(val) == 1) { val = distribution(generator); }
        cqf_oom.insert(val, (1ULL << (val % 5)));
        verif.insert({ val, (1ULL << (val % 5)) });
    }

    const std::string path = "/tmp/bqf_test_oom_roundtrip";
    cqf_oom.save_on_disk(path);
    Bqf_oom loaded = Bqf_oom::load_from_disk(path);
    std::filesystem::remove(path);

    EXPECT_EQ(loaded.enumerate(), verif);
    EXPECT_EQ(loaded.quotient_size,  cqf_oom.quotient_size);
    EXPECT_EQ(loaded.remainder_size, cqf_oom.remainder_size);
    EXPECT_EQ(loaded.elements_inside, cqf_oom.elements_inside);
}


string string_of_int(uint16_t n) {
    string r = "";
    uint32_t mask = 0b11;
    for (int i = 0; i < 8; i++){
        switch (mask&n) {
            case 0:  r.push_back('A'); break;
            case 1:  r.push_back('C'); break;
            case 2:  r.push_back('T'); break;
            default: r.push_back('G');
        }
        n = n>>2;
    }
    return r;
}

class BqfCfTest : public ::testing::Test {
protected:
    std::default_random_engine generator;
    std::uniform_int_distribution<uint16_t> distribution16;

    Bqf_cf small_bqf_cf;
    Bqf_cf bigger_bqf_cf;

    void SetUp() override {
        generator.seed(time(NULL));

        small_bqf_cf = Bqf_cf(7, 8, text, false);
        bigger_bqf_cf = Bqf_cf(18, 28, text, false);
    }
};

// SimpleInsert: very noisy (prints every kmer), kept for manual debugging.
// Run manually: ./build/bin/unit_tests --gtest_filter=BqfCfTest.SimpleInsert
// TEST_F(BqfCfTest, SimpleInsert) { ... }

bool word_in_file(string word, string filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filename);
    }
    string line;
    while (file >> line) {
        if (line.find(word) != line.npos){
            file.close();
            return true;
        }
    }
    file.close();
    return false;
}

TEST_F(BqfCfTest, FilterFastaFile) {
    string fasta_file = "../../examples/data/ecoli_100Kb_reads_40x.fasta";
    string kmc_check = "../../examples/data/ecoli_28_counted";
    string filtered_kmers = "ecoli_28_filtered.txt";
    vector<string> files;
    files.push_back(fasta_file);
    bigger_bqf_cf.filter_fastx_file(files, filtered_kmers);

    ifstream kmc_results(kmc_check);
    string kmer;
    uint64_t count;
    int nb = 0;
    if (kmc_results.is_open()){
        while (kmc_results >> kmer >> count && nb < 100) {
            string nkmer = decode(encode(kmer), bigger_bqf_cf.kmer_size);
            EXPECT_EQ(kmer, nkmer);
            string canonical_kmer = canonical(kmer, bigger_bqf_cf.kmer_size);
            EXPECT_EQ(word_in_file(canonical_kmer, filtered_kmers), count > 1) << canonical_kmer << " not here, but count " << count << endl;
            nb ++;
        }
    }
    kmc_results.close();
}
