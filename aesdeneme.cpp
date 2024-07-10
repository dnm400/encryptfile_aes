#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <cstdint> // For uint8_t
#include <algorithm>
#include <cstring>
#include <cctype> // for isspace


#define Nr 10 //number of rounds
#define Nk 4 //key length
#define Nb 4 //block size

using namespace std;
static const uint8_t sbox[256] = {
  //0     1    2      3     4    5     6     7      8    9     A      B    C     D     E     F
  0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76, //0
  0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0, //1
  0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15, //2 
  0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75, //3
  0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84, //4
  0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf, //5
  0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8, //6
  0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2, //7
  0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73, //8
  0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb, //9
  0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79, //A
  0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08, //B
  0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a, //C
  0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e, //D
  0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf, //E
  0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 }; //F

static const uint8_t inversesbox[256] = {
  0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,
  0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,
  0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
  0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,
  0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,
  0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
  0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,
  0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,
  0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
  0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,
  0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,
  0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
  0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,
  0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,
  0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
  0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d };

void SubBytes(vector<vector<uint8_t>>& subst){ 
for (int i = 0; i < 4; ++i){
    for(int j = 0; j < 4; ++j){
        subst[i][j] = sbox[subst[i][j]];
    }
    }
 }


void invSubBytes(vector<vector<uint8_t>>& subst){ 
for (int i = 0; i < 4; ++i){
    for(int j = 0; j < 4; ++j){
        subst[i][j] = inversesbox[subst[i][j]];
    }
    } 
 }
void ShiftRows(vector<vector<uint8_t>>& newst) { //vector of vectors provides matrix
    //second row, one left shift
    vector<vector<uint8_t>> temp = newst;
    newst[0][1] = temp[1][1];
    newst[1][1] = temp[2][1];
    newst[2][1] = temp[3][1];
    newst[3][1] = temp[0][1];
    //third row, two shift
    newst[0][2] = temp[2][2];
    newst[1][2] = temp[3][2];
    newst[2][2] = temp[0][2];
    newst[3][2] = temp[1][2];
    //fourth row, three shift
    newst[0][3] = temp[3][3];
    newst[1][3] = temp[0][3];
    newst[2][3] = temp[1][3];
    newst[3][3] = temp[2][3];
}


void invShiftRows(vector<vector<uint8_t>>& newst) { //vector of vectors provides matrix
    //second row, one right shift
    vector<vector<uint8_t>> temp = newst;
    newst[1][0] = temp[1][3];
    newst[1][1] = temp[1][0];
    newst[1][2] = temp[1][1];
    newst[1][3] = temp[1][2];
    //third row, two shift
    newst[2][0] = temp[2][2];
    newst[2][1] = temp[2][3];
    newst[2][2] = temp[2][0];
    newst[2][3] = temp[2][1];
    //fourth row, three shift
    newst[3][0] = temp[3][1];
    newst[3][1] = temp[3][2];
    newst[3][2] = temp[3][3];
    newst[3][3] = temp[3][0];
}


void RotWord(vector<uint8_t>& rotw){
    uint8_t temp = rotw[0];
    rotw[0] = rotw[1];
    rotw[1] = rotw[2];
    rotw[2] = rotw[3];
    rotw[3] = temp;
}

uint8_t GalF(uint8_t var1, uint8_t var2) { //ChatGPT Version
    uint8_t result = 0x00;
    for (int i = 0; i < 8; ++i) {
        if (var2 & 1) {
            result ^= var1;
        }
        bool hi_bit_set = (var1 & 0x80);
        var1 <<= 1;
        if (hi_bit_set) {
            var1 ^= 0x1b; // XOR with the irreducible polynomial
        }
        var2 >>= 1;
    }
    return result;
}


void MixColumns(vector<uint8_t>& mixc){ 
    //get vectors of matrix
    vector<uint8_t> temp = mixc; 
    mixc[0] = GalF(0x02, temp[0]) ^ GalF(0x03, temp[1]) ^ temp[2] ^ temp[3];
    mixc[1] = temp[0] ^ GalF(0x02, temp[1]) ^ GalF(0x03, temp[2]) ^ temp[3];
    mixc[2] = temp[0] ^ temp[1] ^ GalF(0x02, temp[2]) ^ GalF(0x03, temp[3]);
    mixc[3] = GalF(0x03, temp[0]) ^ temp[1] ^ temp[2] ^ GalF(0x02, temp[3]);

}

void invMixColumns(vector<uint8_t>& invmix){ 
    //get vectors of matrix
    vector<uint8_t> temp = invmix; 
    invmix[0] = GalF(0x0e, temp[0]) ^ GalF(0x0b, temp[1]) ^  GalF(0x0d, temp[2]) ^  GalF(0x09, temp[3]);
    invmix[1] =  GalF(0x09, temp[0]) ^ GalF(0x0e, temp[1]) ^ GalF(0x0b, temp[2]) ^  GalF(0x0d, temp[3]);
    invmix[2] =  GalF(0x0d, temp[0]) ^  GalF(0x09, temp[1]) ^ GalF(0x0e, temp[2]) ^ GalF(0x0b, temp[3]);
    invmix[3] = GalF(0x0b, temp[0]) ^  GalF(0x0d, temp[1]) ^  GalF(0x09, temp[2]) ^ GalF(0x0e, temp[3]);
}

//For Key Schedule
static const uint8_t Rcon[10] = { //animation video
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36
};

//just for 1 ROUND, assuming 4x4 cipher is given
void updateCipher(vector<vector<uint8_t>>& updatedkey, uint8_t Rcont ){
    vector<vector<uint8_t>> keyschedule = updatedkey;
    vector<uint8_t> veccolumn(4);
    for(int i = 0; i < 4; ++i){ 
        veccolumn[i] = keyschedule[3][i];
    }
    RotWord(veccolumn);
    for(int i = 0; i < 4; ++i){
        veccolumn[i] = sbox[veccolumn[i]];
        }
        veccolumn[0] ^= Rcont;
    for(int i = 0; i < 4; ++i){   
        updatedkey[0][i] = keyschedule[0][i] ^ veccolumn[i];
        updatedkey[1][i] = updatedkey[0][i] ^ keyschedule[1][i];
        updatedkey[2][i] = updatedkey[1][i] ^ keyschedule[2][i];
        updatedkey[3][i] = updatedkey[2][i] ^ keyschedule[3][i];
    }

}  

void AddRoundKey(vector<vector<uint8_t>>& state, vector<vector<uint8_t>>& RoundKey ){ //same function for inverse version
    for(int i = 0; i < 4; ++i){ 
        for(int j = 0; j < 4; ++j){
            state[i][j] ^= RoundKey[i][j];
        }

    } 
} 


void strtomat(const string& hexString, vector<vector<uint8_t>>& matrix) {

    // Convert each pair of hex digits and fill the matrix column by column
    for (int i = 0; i < 16; ++i) {
        char hexPair[3] = { hexString[2*i], hexString[2*i + 1], '\0' };
        matrix[i / 4][i % 4] = static_cast<uint8_t>(strtoul(hexPair, nullptr, 16));
    }
}


string mattostr(vector<vector<uint8_t>> &mat) {
    stringstream ss;
    for (const auto &row : mat) {
        for (const auto &val : row) {
            ss << hex << setw(2) << setfill('0') << static_cast<int>(val);
        }
    }
    return ss.str();
}

void CTRtomat(uint8_t CTR[16], vector<vector<uint8_t>> &CTRmat) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            CTRmat[j][i] = CTR[j * 4 + i];
        }
    }
}

void incrementCTR(uint8_t CTR[16]) {
    for (int i = 15; i >= 12; --i) { // Only increment the counter part (last 4 bytes)
        if (CTR[i] != 0xff) {
            CTR[i] += 1;
            break;
        } else {
            CTR[i] = 0x00; 
            CTR[i-1] += 1;
        }
    }
}


int main(){ //define types

    uint8_t IV[12]= {0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xaa, 0xbb}; // 96-bit IV
    uint8_t counter32[4] = {0x00, 0x00, 0x00, 0x00}; //remaining 32-bit
    uint8_t CTR[sizeof(IV) + sizeof(counter32)]; //16 byte
    memcpy(CTR, IV, sizeof(IV));
    memcpy(CTR + sizeof(IV), counter32, sizeof(counter32));
    vector<vector<uint8_t>> key(4, vector<uint8_t>(4));

    cout << "Plain Text? " <<endl;
    string plainin;
    getline(cin, plainin);
    plainin.erase(remove_if(plainin.begin(), plainin.end(), [](char c) { return isspace(static_cast<unsigned char>(c)); }), plainin.end());

    size_t numBlocks = plainin.length() / 32; //blocks for plaintext


    cout << "Key? " << endl;
    string keyin;
    getline(cin, keyin);
    keyin.erase(remove_if(keyin.begin(), keyin.end(), [](char c) { return isspace(static_cast<unsigned char>(c)); }), keyin.end());


 for (size_t m = 0; m < numBlocks; ++m) {
    string block = plainin.substr(m * 32, 32); // plaintext shouldnt be written seperate
    vector<vector<uint8_t>> plaintext(4, vector<uint8_t>(4));
    strtomat(block, plaintext);
    vector<vector<uint8_t>> CTRmat(4, vector<uint8_t>(4)); //matrix of counter
    CTRtomat(CTR, CTRmat);
    strtomat(keyin, key);
    AddRoundKey(CTRmat,key);


    for(int i = 1; i < Nr; ++i){ //Number of rounds Nr
        SubBytes(CTRmat);
        ShiftRows(CTRmat);
        vector<uint8_t> vec0(4), vec1(4), vec2(4), vec3(4);
        for(int j = 0; j < 4; ++j){ 
            vec0[j] = CTRmat[0][j];
            vec1[j] = CTRmat[1][j];
            vec2[j] = CTRmat[2][j];
            vec3[j] = CTRmat[3][j];
        }
        MixColumns(vec0);
        MixColumns(vec1);
        MixColumns(vec2);
        MixColumns(vec3);

        for(int j = 0; j < 4; ++j){ 
            CTRmat[0][j] = vec0[j] ;
            CTRmat[1][j] = vec1[j];
            CTRmat[2][j] = vec2[j];
            CTRmat[3][j] = vec3[j];
        }

        updateCipher(key, Rcon[i - 1]);
        AddRoundKey(CTRmat, key);
    }   

    updateCipher(key, Rcon[Nr-1]);

    //for last round no mixcolumn 
    SubBytes(CTRmat);
    ShiftRows(CTRmat);
    AddRoundKey(CTRmat, key);

    for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                plaintext[j][k] ^= CTRmat[j][k];
            }
        }
   
    string ciphertext = mattostr(plaintext);
    cout << ciphertext ;

    incrementCTR(CTR);
    }
    cout << " is the ciphertext.";
    return 0;

}