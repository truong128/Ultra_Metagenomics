# Ultra_Metagenomics

The public guide is `README.md`. This file keeps the earlier Vietnamese operator notes.

# MangroveMetaX Ultra v3.1

**Default profile is now `linux64` for this workstation (Ryzen 9 9950X3D, 64 GB RAM, NVIDIA RTX PRO 5000 72 GB).** One command still runs every computational stage, records BLOCKED tools/databases, and keeps duplicate engines out of the production tables. The Mac 64 GB profile is preserved in `config/project_mac64.json`.

```bash
cd MetagenomicsX
python3 metax.py doctor
python3 metax.py run-all --dry-run
python3 metax.py run-all
```

Read `docs/LINUX64_SETUP.md` for the 48 GB container cap and the swap-backed GTDB-Tk stage, and `docs/PIPELINE_INTEGRATION.md` for how SqueezeMeta, nf-core/taxprofiler, DRAM, METABOLIC, HUMAnN, metagenome-atlas and the AAFC host-associated workflow were merged without double-counting.

# MangroveMetaX Ultra v3.0 AutoRun

**Một lệnh để điều phối tất cả các adapter phân tích tin sinh học đã triển khai, ghi rõ các nhánh còn thiếu.** Bản Mac 64 GB bên dưới vẫn đúng khi `resource_mode` là `mac64`. **Không yêu cầu GPU ở chế độ Mac.**

## Chạy một lần (sau khi cài các công cụ, database và nhập FASTQ thật)

```bash
cd MangroveMetaX_Ultra_v3.0_autorun
# Nhập R1/R2 thực tế và metadata vào config/samples.tsv
# Cấu hình reference_databases + core flags trong config/project.json
python3 metax.py doctor
python3 metax.py run-all --dry-run
python3 metax.py run-all
```

Chạy lại cùng lệnh sau khi mất kết nối hoặc dừng: những stage có dấu hiệu output hợp lệ được tái sử dụng; nf-core/mag cũng có `-resume`. Để buộc chạy mới toàn bộ stage: `python3 metax.py run-all --force`. Xem báo cáo cuối tại `results/execution/report.html`, `report.json`, `status.json` và logs từng stage.

**QUAN TRỌNG:** Lệnh `run-all` mặc định trả về **mã lỗi** khi bất kỳ bước nào thất bại, bị chặn vì thiếu công cụ/database, hoặc không thể thực hiện do không có MAG đủ QC. Nếu muốn *cố ý* tạo báo cáo từ phần kết quả khả dụng, dùng `--allow-partial`, nhưng đây không phải kết quả full end-to-end.

## Các nhánh tự điều phối

1. FASTQ integrity (toàn bộ reads).
2. nf-core/mag v5.5.0, QC, MEGAHIT assembly, Prodigal, mapping, MetaBAT2/MaxBin2, DAS Tool, CheckM2; nếu kích hoạt đúng flags/database thì thêm geNomad, MetaEuk, GUNC.
3. Bridge FASTA, ORF, MAG candidates và audit assembly.
4. Auto MAG screening theo CheckM2 >=50% completeness và <=10% contamination (khi match ID chính xác). Không tự nhận định MAG đã được chuyên gia duyệt hay nhận diện loài.
5. skani + GTDB representative sketches (chỉ ANI candidates; chưa phải GTDB-Tk species assignment).
6. Kraken2, Bracken, phyloFlash, Tiara.
7. eggNOG-mapper, dbCAN, DRAM2, FeGenie, AMRFinderPlus, antiSMASH, METABOLIC-C.
8. Tổng hợp bảng KO từ eggNOG (gene counts, chưa phải pathway activity), ba hình SVG mô tả assembly từ FASTA.
9. CheckV sau geNomad; báo cáo kết quả thật, provenance và các bước bị chặn.

Đây là **21 computational execution stages** bao phủ các tác vụ lõi và nhánh triển khai được; bảng nghiên cứu 91 bước tại `docs/STEP_REGISTRY.tsv` còn bao gồm các bước **thiết kế nghiên cứu, lấy mẫu, wet-lab, kiểm chứng phylogeny chuyên sâu, thống kê cần nhiều mẫu và đánh giá của chuyên gia**. Những bước đó không thể được xác nhận hoàn thành bằng một lệnh khi chưa có dữ liệu hoặc công việc thực nghiệm. Báo cáo tự ghi sự khác biệt này.

## Vì sao có thể xuất hiện BLOCKED?

Các công cụ và reference databases lớn **không được cung cấp trong ZIP**. Chúng cần được cài đặt riêng, và phải khai báo đường dẫn thật trong `config/project.json`. Không tự tải database hàng trăm GB mà không biết dung lượng SSD và điều khoản phân phối. Các bước nặng như GTDB-Tk pplacer vẫn dành cho `hpc`, **không được khởi chạy trong `mac64`**. Chương trình không bịa kết quả nếu một công cụ thiếu hoặc sai phiên bản.

Để các nhánh virus và nhân thực của nf-core/mag chạy, cấu hình `core_enable_genomad: true` với `reference_databases.genomad`, và `core_enable_metaeuk: true` với `reference_databases.metaeuk`. Chạy `python3 metax.py core --dry-run` để kiểm tra lệnh chính trước khi thực thi. Nếu không có viral contigs, CheckV có thể không có đối tượng để phân tích. Cung cấp `metabolic_c_script` nếu dùng METABOLIC-C và lưu ý bước này cần chỗ trống để giải nén FASTQ.

Mặc định `config/samples.tsv` **chứa đường dẫn minh họa**, không phải dữ liệu thật. ZIP này **không chứa FASTQ**, không cài database, không thực hiện phép tính sinh học trên dữ liệu Cần Giờ. Các kiểm thử local kiểm tra logic, không xác nhận toàn bộ toolchains trên Apple Silicon.

## Tài liệu

- `docs/AUTORUN.md`: command, status codes, database, QC, interpretation.
- `docs/MAC64_SETUP.md`: Docker 48 GB, Nextflow 40 GB/task, Apple Silicon caveats.
- `docs/STEP_REGISTRY.tsv`: toàn bộ 91 bước nghiên cứu, kể cả bước cần kiểm chứng ngoài pipeline.
- `docs/QUALITY_GATES.md`: các ngưỡng QC và yêu cầu bằng chứng.
- `docs/MARKER_PANEL.tsv`: hệ thống diagnostic functional markers.
- `docs/DATABASE_SETUP.md`: nguồn và điều kiện cài đặt databases.

## Unit tests

```bash
python3 -m unittest discover -s tests -v
```

**Tài liệu kỹ thuật gốc:** https://nf-co.re/mag/5.5.0/ ; https://github.com/bluenote-1577/skani ; https://github.com/BortonWrightonLabs/DRAM ; https://github.com/AnantharamanLab/METABOLIC ; https://run-dbcan.readthedocs.io/

## Citation

Truong Nguyen, CUNY School of Medicine.
