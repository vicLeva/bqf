#include <gtest/gtest.h>

#include "filter.hpp"
#include "bcqf_ec.hpp" 
#include "bcqf_oom.hpp"
#include <random>

using namespace std;

class RsqfTest : public ::testing::Test {
 protected:
  void SetUp() override {
    generator.seed(time(NULL));
    
    usual_qf = Rsqf(2);
    small_qf = Rsqf(7, 64-7, false);
    cqf = Bcqf_ec(4, 5, false);
  }

  // void TearDown() override {}

  std::default_random_engine generator;
  std::uniform_int_distribution<uint64_t> distribution;

  Rsqf usual_qf;
  Rsqf small_qf;
  Bcqf_ec cqf;
};


/* TEST_F(RsqfTest, QuotientSize) {
    EXPECT_EQ(usual_qf.get_quot_size(), 19);
}


TEST_F(RsqfTest, get_runstart) {
    small_qf.insert((2ULL<<30)+64);
    EXPECT_EQ(small_qf.get_runstart(64), 64);
    small_qf.insert((2ULL<<31)+63);

    small_qf.insert((2ULL<<32)+126);
    small_qf.insert((2ULL<<33)+126);
    small_qf.insert((2ULL<<34)+126);
    EXPECT_EQ(small_qf.get_runstart(64), 64);
    std::cout << small_qf.block2string(0) << "\n" << small_qf.block2string(1);
}



TEST_F(RsqfTest, insert_pos_0) {
    small_qf.insert((2ULL << 30) + 90); 
    small_qf.insert(2ULL << 31); 
    //std::cout << small_qf.block2string(0) << "\n" << small_qf.block2string(1);
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

TEST_F(RsqfTest, offset2_runOverflowsBy1) {
    small_qf.insert((3ULL<<31) + 64);

    small_qf.insert((3ULL<<30) + 63);
    small_qf.insert((3ULL<<30) + 63);


    EXPECT_EQ(small_qf.get_remainder(64), 25165824ULL);
}

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
    //run of quot 20, starts @20 ends up @79
    for (int i = 0; i < 69; i++){ usual_qf.insert((1ULL<<32)+ 40); } 
    //run of quot 40, starts @80 ends up @148
    usual_qf.insert((1ULL<<32)+ 149);
    //run of quot 149, starts @149 ends up @149 
    //(boundary shift deletion in case we delete in 2nd run, because can't be shifted)
    for (int i = 0; i < 11; i++){ usual_qf.insert((1ULL<<32)+ 150); }
    //run of quot 150, starts @150 ends up @160 (FUS == 161)

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
    //bugs because of offset2_runOverflowsBy1
    std::unordered_set<uint64_t> verif;
    
    small_qf.insert((3ULL<<30) + 35); verif.insert((3ULL<<30) + 35);
    small_qf.insert((3ULL<<31) + 35); verif.insert((3ULL<<31) + 35);
    small_qf.insert((3ULL<<25) + 35); verif.insert((3ULL<<25) + 35);

    small_qf.insert((3ULL<<31) + 64); verif.insert((3ULL<<31) + 64);

    small_qf.insert((3ULL<<30) + 63); verif.insert((3ULL<<30) + 63);
    small_qf.insert((3ULL<<30) + 63); verif.insert((3ULL<<30) + 63);


    EXPECT_EQ(small_qf.enumerate(), verif);
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
    //run of quot 20, starts @20 ends up @35
    for (int i = 0; i < 28; i++){ small_qf.insert((40ULL<<7)+ 40); } 
    //run of quot 40, starts @40 ends up @67
    small_qf.insert((99ULL<<7)+ 99);
    //run of quot 99, starts @99 ends up @99 
    for (int i = 0; i < 12; i++){ small_qf.insert((100ULL<<7)+ 100); }
    //run of quot 100, starts @100 ends up @111 
    
    for (int i = 0; i < 48; i++){ small_qf.insert((96ULL<<7)+ 96); }
    //run of quot 96, starts @96 ends up @15
    //run of quot 99 is shifted, starts @16 ends up @16 
    //run of quot 100 is shifted, starts @17 ends up @28 
    //run of quot 20 is shifted, starts @29 ends up @44
    //run of quot 40 is shifted, starts @45 ends up @72


    compare = std::make_pair(29, 44);   EXPECT_EQ(small_qf.get_run_boundaries(20), compare);
    compare = std::make_pair(45, 72);   EXPECT_EQ(small_qf.get_run_boundaries(40), compare);
    compare = std::make_pair(96, 15);   EXPECT_EQ(small_qf.get_run_boundaries(96), compare);
    compare = std::make_pair(16, 16);   EXPECT_EQ(small_qf.get_run_boundaries(99), compare);
    compare = std::make_pair(17, 28);   EXPECT_EQ(small_qf.get_run_boundaries(100), compare);

    
    std::cout << small_qf.block2string(0) << "\n" << small_qf.block2string(1);


    for (int i = 0; i < 28; i++){ small_qf.remove((96ULL<<7)+ 96); }

    std::cout << small_qf.block2string(0) << "\n" << small_qf.block2string(1);


    compare = std::make_pair(20, 35);   EXPECT_EQ(small_qf.get_run_boundaries(20), compare);
    compare = std::make_pair(40, 67);   EXPECT_EQ(small_qf.get_run_boundaries(40), compare);
    compare = std::make_pair(96, 115);   EXPECT_EQ(small_qf.get_run_boundaries(96), compare);
    compare = std::make_pair(116, 116);   EXPECT_EQ(small_qf.get_run_boundaries(99), compare);
    compare = std::make_pair(117, 0);   EXPECT_EQ(small_qf.get_run_boundaries(100), compare);
}



TEST_F(RsqfTest, first_unused_slot1) {
    for (int i = 0; i < 2; i++){ small_qf.insert((123<<7) + 123); }
    small_qf.insert((124<<7) + 124);
    for (int i = 0; i < 2; i++){ small_qf.insert((125<<7) + 125); }
    small_qf.insert((126<<7) + 126);
    small_qf.insert((1<<7) + 1);
    for (int i = 0; i < 3; i++){ small_qf.insert((2<<7) + 2); }

    std::cout << small_qf.block2string(0) << "\n" << small_qf.block2string(1);


    EXPECT_EQ(small_qf.first_unused_slot(123), 5);
}

TEST_F(RsqfTest, first_unused_slot2) {
    for (int i = 0; i < 2; i++){ small_qf.insert((123<<7) + 123); }
    for (int i = 0; i < 4; i++){ small_qf.insert((124<<7) + 124); }
    for (int i = 0; i < 3; i++){ small_qf.insert((1<<7) + 1); }

    std::cout << small_qf.block2string(0) << "\n" << small_qf.block2string(1);


    EXPECT_EQ(small_qf.first_unused_slot(123), 4);
}


TEST_F(RsqfTest, first_unused_slot3) {
    small_qf.insert((63<<7) + 63);
    small_qf.insert((65<<7) + 65);

    std::cout << small_qf.block2string(0) << "\n" << small_qf.block2string(1);


    EXPECT_EQ(small_qf.first_unused_slot(63), 64);
}


TEST_F(RsqfTest, shift_bits_right_metadata) {
    small_qf.insert((63<<7) + 63);
    small_qf.insert((3000<<7) + 63);
    small_qf.insert((64<<7) + 64);

    small_qf.remove((3000<<7) + 63);

    std::cout << small_qf.block2string(0) << "\n" << small_qf.block2string(1);

    EXPECT_EQ(small_qf.get_runend_word(0), 1ULL<<63);
} */

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

TEST_F(RsqfTest, finalTest) {
    uint64_t val;
    unordered_set<uint64_t> verif; 

    //INSERT
    for (uint64_t i=0; i < (1ULL<<18)-1; i++){
        val = distribution(generator);      
        usual_qf.insert(val);
        verif.insert(val);
    }

    EXPECT_EQ(usual_qf.enumerate(), verif);

    //REMOVE
    for (uint64_t i=0; i < (1ULL<<18)-1; i++){
        val = *verif.begin(); 
        verif.extract(val);    
        usual_qf.remove(val);
    }

    EXPECT_EQ(usual_qf.enumerate(), verif);
}






class BCqfTest : public ::testing::Test {
 protected:
  void SetUp() override {
    generator.seed(time(NULL));
    
    small_cqf = Bcqf_ec(7, 64-7, 5, false);
    cqf = Bcqf_ec(1, 5, false);
    
    small_cqf_oom = Bcqf_oom(7, 64-7, 5, false);
    cqf_oom = Bcqf_oom(1, 5, false);
  }

  // void TearDown() override {}

  std::default_random_engine generator;
  std::uniform_int_distribution<uint64_t> distribution;

  Bcqf_ec small_cqf;
  Bcqf_ec cqf;
  
  Bcqf_oom small_cqf_oom;
  Bcqf_oom cqf_oom;
};


TEST_F(BCqfTest, insert1occs) {
    uint64_t val;
    std::map<uint64_t, uint64_t> verif; 

    //INSERT
    for (uint64_t i=0; i < (1ULL<<17)-1; i++){
        val = distribution(generator);    
        while (verif.count(val) == 1) { //already seen key
            val = distribution(generator);    
        }  

        cqf.insert(val);
        verif.insert({ val, 1 });
    }

    EXPECT_EQ(cqf.enumerate(), verif);

    //REMOVE
    std::map<uint64_t,uint64_t>::iterator it;
    for (it = verif.begin(); it != verif.end(); it++){
        cqf.remove((*it).first);
    }
    verif.clear();

    EXPECT_EQ(cqf.enumerate(), verif);
}


TEST_F(BCqfTest, insertRDMoccs) {
    uint64_t val;
    std::map<uint64_t, uint64_t> verif; 

    //INSERT
    for (uint64_t i=0; i < (1ULL<<17)-1; i++){
        val = distribution(generator);    
        while (verif.count(val) == 1) { //already seen key
            val = distribution(generator);    
        }  

        cqf.insert(val, val%31);
        verif.insert({ val, val%31 });
    }

    EXPECT_EQ(cqf.enumerate(), verif);

    //REMOVE
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

    //INSERT
    for (uint64_t i=0; i < (1ULL<<17)-1; i++){
        val = distribution(generator);    
        while (verif.count(val) == 1) { //already seen key (unlikely)
            val = distribution(generator);    
        }  
        
        cqf_oom.insert(val, (1ULL << val%31));
        verif.insert({ val, (1ULL << val%31) });
    }

    EXPECT_EQ(cqf_oom.enumerate(), verif);

    //REMOVE
    std::map<uint64_t,uint64_t>::iterator it;
    for (it = verif.begin(); it != verif.end(); it++){
        cqf_oom.remove((*it).first);
    }
    verif.clear();

    EXPECT_EQ(cqf_oom.enumerate(), verif);
}