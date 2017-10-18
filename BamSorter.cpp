#include <array>
#include <cstdio>
#include <cstdint>
#include <future>
#include <iostream>
#include <sstream>

#include "boost/filesystem.hpp"
#include "boost/variant.hpp"


#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <lz4.h>
#include <zlib.h>

#include <string_view>

struct DecompressedData {
    DecompressedData(boost::variant<std::string, std::string_view> decompressed_data) : decompressed_data_(std::move(decompressed_data)) {}

    boost::variant<std::string, std::string_view> decompressed_data_;
};

struct CompressedBlock {
    CompressedBlock(boost::variant<std::string, std::string_view> compressed_block) : compressed_block_(std::move(compressed_block)) {}

    boost::variant<std::string, std::string_view> compressed_block_;
};

struct string_view_visitor : public boost::static_visitor<std::string_view> {
    std::string_view operator() (std::string_view s) const { return s; }
    std::string_view operator() (const std::string & s) const { return s; }
};

class BgzfBlock {
public:
    BgzfBlock(CompressedBlock compressed_block) : decompressed_data_(std::string()), compressed_block_(std::move(compressed_block)) {}

    std::string decompress_zlib(std::string_view compressed_data, uint32_t isize) {
        std::string decompressed_data;
        decompressed_data.resize(isize);

        z_stream zs;
        zs.next_in = (Bytef *)(compressed_data.begin());
        zs.avail_in  = compressed_data.size();
        zs.next_out  =  (Bytef *)decompressed_data.data();
        zs.avail_out = decompressed_data.size();
        zs.zalloc = nullptr;
        zs.zfree = nullptr;

        int status = inflateInit2(&zs, -15);
        if (status != Z_OK) {
            throw std::runtime_error("Error initializing zlib");
        }

        status = inflate(&zs, Z_FINISH);
        if ( status != Z_STREAM_END ) {
            inflateEnd(&zs);
            throw std::runtime_error("Zlib inflate failed");
        }

        assert(zs.total_out == isize);

        status = inflateEnd(&zs);
        if (status != Z_OK) {
            throw std::runtime_error("ending zlib infiate failed");
        }

        return decompressed_data;
    }

    std::string decompress_lz4(std::string_view compressed_data, uint32_t isize) {
        std::string decompressed_data;
        decompressed_data.resize(isize);
        int decompressed_size = LZ4_decompress_safe(compressed_data.data(), decompressed_data.data(), compressed_data.size(), isize);
        if (decompressed_size < 0) {
            throw std::runtime_error("lz4 failed to decompress");
        }

        if (decompressed_size != isize) {
            throw std::runtime_error("Decompressed size does not match what we expected from reading the BGZF header?");
        }

        return decompressed_data;
    }

    uint16_t compress_zlib(const std::string_view & decompressed_data, std::string & compressed_block) {
        z_stream zs;
        zs.next_in   = (Bytef*)decompressed_data.begin();
        zs.avail_in  = decompressed_data.size();
        zs.next_out  = (Bytef*)&compressed_block[18];
        zs.avail_out = 65536 - 18 - 8;
        zs.zalloc = nullptr;
        zs.zfree = nullptr;
        int status = deflateInit2(&zs, 6, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY);
        if (status != Z_OK) {
            throw std::runtime_error("Failed to initialize Zlib");
        }

        status = deflate(&zs, Z_FINISH);
        if (status != Z_STREAM_END) {
            throw std::runtime_error("Failed to deflate data");
        }

        assert(zs.total_out <= 65536 - 18 - 8);

        return zs.total_out + 18 + 8 - 1;
    }

    uint16_t compress_lz4(const std::string_view & decompressed_data, std::string & compressed_block) {
        int compressed_size = LZ4_compress_default(decompressed_data.data(), compressed_block.data() + 18, decompressed_data.size(), 65536-18-8);
        if (compressed_size == 0) {
            throw std::runtime_error("Input (decompressed_data) too large to compress into 65510 bytes");
        }

        if (compressed_size < 0) {
            throw std::runtime_error("Error compressing file");
        }

        assert(compressed_size <= 65536-18-8);

        return compressed_size + 18 + 8 - 1;
    }

    void decompress() {
        const std::string_view & compressed_block = boost::apply_visitor(string_view_visitor(), compressed_block_.compressed_block_);

        auto isize_raw_it = compressed_block.end() - 4;

        uint32_t isize;
        memcpy(&isize, isize_raw_it, sizeof(isize));

        std::string_view compressed_data(compressed_block.begin() + 18, compressed_block.size() - 1 - 6 - 19);

        std::string decompressed_data = decompress_lz4(compressed_data, isize);

        decompressed_data_.decompressed_data_ = std::move(decompressed_data);
    }

    void compress() {
        const std::string_view & decompressed_data = boost::apply_visitor(string_view_visitor(), decompressed_data_.decompressed_data_);

        std::string compressed_block;
        compressed_block.resize(65536);

        compressed_block[0] = 31;   // ID1
        compressed_block[1] = 139;  // ID2
        compressed_block[2] = 8;    // CM
        compressed_block[3] = 4;    // FLG

        // MTIME
        compressed_block[4] = 0;
        compressed_block[5] = 0;
        compressed_block[6] = 0;
        compressed_block[7] = 0;

        compressed_block[8] = 0;    // XFL
        compressed_block[9] = 255;  // OS (255 = unknown OS)

        // XLEN
        compressed_block[10] = 0;
        compressed_block[11] = 0;

        compressed_block[12] = 66;   // SI1
        compressed_block[13] = 67;   // SI2

        // SLEN
        compressed_block[14] = 2;
        compressed_block[15] = 0;

        // CDATA
        uint16_t bsize = compress_lz4(decompressed_data, compressed_block);

        // BSIZE
        memcpy(&compressed_block[16], &bsize, sizeof(bsize));

        // CRC32
        uint32_t crc = crc32(crc32(0, NULL, 0), (Bytef*)&decompressed_data[0], decompressed_data.size());
        memcpy(&compressed_block[bsize+1-8], &crc, sizeof(crc));

        // ISIZE
        uint32_t isize = decompressed_data.size();
        memcpy(&compressed_block[bsize+1-4], &isize, sizeof(isize));

        compressed_block.resize(bsize+1);

        compressed_block_.compressed_block_ = std::move(compressed_block);
    }

    DecompressedData decompressed_data_;
    CompressedBlock compressed_block_;
};

uint64_t get_file_size(int fd) {
    struct stat sb;
    auto res = fstat(fd, &sb);
    if (res != 0) {
        std::stringstream ss;
        ss << "Failed to fstat, fd: " << fd << " " << strerror(errno);
        throw std::runtime_error(ss.str());
    }

    return sb.st_size;
}

uint16_t get_bgzf_block_end_coords(uint64_t start_coord, char * input_file_ptr) {
    uint16_t bsize;
    memcpy(&bsize, input_file_ptr + start_coord + 16, sizeof(bsize));

    return bsize + 1;
}

void write_raw(const boost::filesystem::path & output_path, const std::vector<BgzfBlock> & blocks) {
    std::ofstream output_file(output_path.string(), std::ios::out | std::ios::binary);
    
    for (const auto & block : blocks) {
        const std::string_view & decompressed_data = boost::apply_visitor(string_view_visitor(), block.decompressed_data_.decompressed_data_);
        output_file << decompressed_data;
    }
}

void rewrite(const boost::filesystem::path & output_path, std::vector<BgzfBlock> & blocks) {
    int num_threads = 1;
//    int num_threads = std::thread::hardware_concurrency();

    int num_blocks_per_thread = blocks.size() / num_threads;
    int num_remaining_tasks = blocks.size() % num_threads;

    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    {
        auto begin = blocks.begin();
        for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
            auto end = begin + num_blocks_per_thread;

            if (thread_id < num_remaining_tasks) {
                ++end;
            }

            auto threadproc = [begin, end]() {
                for (auto it = begin; it != end; ++it) {
                    it->compress();
                }
            };

            threads.emplace_back(threadproc);

            begin = end;
        }
    }

    for (auto & thread : threads) {
        thread.join();
    }


    std::ofstream output_file(output_path.string(), std::ios::out | std::ios::binary);

    for (auto & block : blocks) {
        const std::string_view & compressed_block = boost::apply_visitor(string_view_visitor(), block.compressed_block_.compressed_block_);

        output_file.write(compressed_block.data(), compressed_block.size());
    }
}

int main(int argc, char * argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " /path/to/input.bam /path/to/output.bam" << std::endl;
        return -1;
    }

    boost::filesystem::path input_file_path(argv[1]);
    boost::filesystem::path output_file_path(argv[2]);

    int input_fd = open(input_file_path.c_str(), O_RDONLY);
    if (input_fd == -1) {
        std::stringstream ss;
        ss << "Failed to open " << input_file_path << " " << strerror(errno);

        throw std::runtime_error(ss.str());
    }

    uint64_t input_file_size = get_file_size(input_fd);
    std::cout << "input file size: " << input_file_size << std::endl;

    void * input_file_mmap_ptr = mmap(NULL, input_file_size, PROT_READ, MAP_PRIVATE, input_fd, 0);
    if (input_file_mmap_ptr == MAP_FAILED) {
        std::stringstream ss;
        ss << "Failed to mmap fd: " << input_fd << " " << strerror(errno);
        throw std::runtime_error(ss.str());
    }

    char * input_file_ptr = reinterpret_cast<char *>(input_file_mmap_ptr);

//    std::vector<char> input_file_vec;
//    input_file_vec.resize(input_file_size);
//    FILE * input_file = fopen(input_file_path.c_str(), "rb");
//    auto num_read = fread(input_file_vec.data(), 1, input_file_size, input_file);
//    if (num_read != input_file_size) {
//        throw std::runtime_error("Number of bytes read not equal to file size");
//    }
//    fclose(input_file);
//
//    char * input_file_ptr = input_file_vec.data();

    uint64_t start_coord = 0;

    std::vector<BgzfBlock> blocks;
    while (true) {
        uint16_t block_size = get_bgzf_block_end_coords(start_coord, input_file_ptr);

        assert(start_coord + block_size - 1 <= input_file_size);

        std::cout << start_coord << " " << block_size << std::endl;

        blocks.emplace_back(CompressedBlock(std::string(input_file_ptr + start_coord, block_size)));

        start_coord += block_size;

        if (start_coord == input_file_size) {
            break;
        }
    }

    std::cout << blocks.size() << " BgzfBlocks detected" << std::endl;

    int num_threads = std::thread::hardware_concurrency();
//    int num_threads = 1;
    int num_trials = 10;

    int num_blocks_per_thread = blocks.size() / num_threads;
    int num_remaining_tasks = blocks.size() % num_threads;

    auto start = std::chrono::high_resolution_clock::now();

    for (int k = 0; k < num_trials; ++k) {
        std::vector<std::thread> threads;
        threads.reserve(num_threads);

        {
            auto begin = blocks.begin();
            for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
                auto end = begin + num_blocks_per_thread;

                if (thread_id < num_remaining_tasks) {
                    ++end;
                }

                auto threadproc = [begin, end]() {
                    for (auto it = begin; it != end; ++it) {
                        it->decompress();
                    }
                };

                threads.emplace_back(threadproc);

                begin = end;
            }
        }

        for (auto & thread : threads) {
            thread.join();
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto decompression_duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    std::cout << "Decompression took: " << decompression_duration_ms << " milliseconds (" << input_file_size * num_trials / (1024 * 1024) / (decompression_duration_ms / 1000.0) << " MiB/s)" << std::endl;

    write_raw(output_file_path, blocks);
    //rewrite(output_file_path, blocks);

    return 0;
}
