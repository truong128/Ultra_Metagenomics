# One-command execution

v3.1 default is `linux64`. See `docs/LINUX64_SETUP.md`. The Mac notes below apply when `resource_mode` is `mac64`.

# One-command execution: v3.0

## Công cụ nền

Python 3.10+, Nextflow, Docker Desktop/Linux Docker, JDK phù hợp với bản Nextflow cài đặt, truy cập mạng để pull containers. External CLI cần có trong PATH của môi trường mà Python điều phối: `skani`, `kraken2`, `bracken`, `phyloFlash.pl`, `tiara`, `emapper.py`, `run_dbcan`, `FeGenie.py`, `amrfinder`, `antismash`, `checkv`, `perl`, cùng package METABOLIC-C. Không mặc định các package đều cài native ARM64. Một số phần phù hợp hơn với Linux x86_64/Conda containers.

**Mac 64 GB:** 8 CPU, Docker 48 GB, max task request 40 GB, một nhiệm vụ lớn một lần. Không gọi GTDB-Tk/pplacer trong `mac64`. Giới hạn yêu cầu tài nguyên không bảo đảm peak RSS thực tế; giám sát Activity Monitor/htop. Không dùng GPU.

## Cấu hình

1. `config/samples.tsv`: sample ID, R1/R2 **đường dẫn tuyệt đối**, giữ barcode/kênh dữ liệu nhất quán, điều kiện mẫu trắng.
2. `config/project.json`: resource_mode=mac64; `reference_databases.*` chỉ định đường dẫn database đang có; `metabolic_c_script` là tệp Perl chính từ bản upstream đã cài; `viral_fasta` để ghi đè khi cần.
3. Kích hoạt `core_enable_genomad` và `core_enable_metaeuk` khi có dữ liệu tham chiếu tương ứng, nếu muốn các công cụ này chạy trong backbone; kiểm tra tham số upstream trước khi chạy phiên bản mới.
4. Không đặt `core_enable_gtdbtk=true` trong mac64. `taxonomy-ani` tự dùng skani, `taxonomy-hpc` chỉ chạy trên máy RAM cao khi được yêu cầu riêng.
5. `python3 metax.py run-all --dry-run` để xem 21 stage và lỗi dependencies. `python3 metax.py run-all` sẽ chạy tuần tự, không tự động bỏ qua.

## Tình trạng

- `SUCCEEDED`: bước thực thi kết thúc và kiểm tra output dự kiến ở mức adapter.
- `CACHED`: cùng cấu hình và metadata FASTQ, đầu ra của stage còn hiện hữu, được dùng lại.
- `FAILED`: lỗi công cụ hoặc lỗi input/output, đọc `results/execution/<stage>.log`.
- `BLOCKED`: thiếu công cụ/database hoặc bước phụ thuộc chưa hoàn thành.
- `NOT_APPLICABLE`: không có MAG hoặc MAG đạt ngưỡng QC theo CheckM2, nên không được suy diễn kết quả ở MAG level.

**Auto MAG gate**: chỉ sử dụng bảng `GenomeBinning/QC/checkm2_summary.tsv`, phải match chính xác MAG ID, completeness >= 50 và contamination <= 10. Fail closed khi không tìm được cột hoặc ID. Không giả định GUNC đã chạy, không tự gán species bằng ngưỡng ANI cố định. Kết quả tự động không thay thế manual validation để xuất bản.

## Giới hạn quan trọng

- DAG này chạy tất cả các adapter tin sinh học được liệt kê, không phải toàn bộ 91 bước nghiên cứu. Các bước nghiên cứu ngoài đời (sample collection, environmental measurements, bench validation, statistics requiring replicates) không thể hoàn thành chỉ bằng FASTQ.
- Protein function ≠ enzyme activity. Co-occurrence ≠ metabolic exchange. MAG ≠ species count. Contig ≠ organism.
- Với một mẫu độc lập, không chạy các kiểm định differential abundance, environmental association và network inference như thể có replicates.
- Source tool versions và database release phải pin/kiểm tra riêng trước khi tuyên bố production validated. ZIP chỉ có unit tests với fixtures nhỏ, chưa có biological integration benchmark.
- Nếu `run-all` kết thúc không thành công, `report.html`/`report.json` vẫn ghi được kết quả nào đã chạy và tại sao các nhánh còn lại không được thực hiện.

## Đọc kết quả

`results/execution/report.html`, `status.json`, `results/audit/audit.json`, `results/summaries/`, `results/figures/`, `results/assets.json`, `results/extensions/*` và `results/taxonomy_ani/ani_candidates.tsv`. `report.html` là **execution report** không phải bản thảo đã được kiểm định khoa học.
