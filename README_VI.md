# Ultra_Metagenomics

Phân tích shotgun metagenome cho mẫu nước. Assembly dùng trong kết quả là nf-core/mag 5.5.0 với MEGAHIT, từng mẫu một. Các bước sau bổ sung phân loại, genome bin và chú giải chức năng. Các bước đó không thay assembly MEGAHIT.

Kho mã gồm workflow và bảng mẫu. Không gồm reads, database tham chiếu, hay kết quả. Reads của ví dụ đã có trên ENA. Bản hướng dẫn đầy đủ bằng tiếng Anh nằm ở `README.md`.

## Ví dụ

Bảng mẫu là hai thư viện hồ Motosu từ nghiên cứu PRJEB35640 (ERP118719). Cả hai được thu vào ngày 7 tháng 9 năm 2015.

| Mẫu | Độ sâu | Run ENA | BioSample |
| --- | --- | --- | --- |
| Motosu_5m | 5 m | ERR14787630 | SAMEA117850005 |
| Motosu_100m | 100 m | ERR14787631 | SAMEA117850006 |

Một ngày và hai độ sâu cho phép so sánh theo độ sâu. Bộ mẫu này không đủ cho một kết luận theo mùa hoặc một mạng đồng xuất hiện.

## Cách chạy

Cài Nextflow, Docker, và những database sẽ dùng. Sao chép profile rồi điền đường dẫn:

```bash
cp config/project.linux64.json config/project.json
python3 metax.py doctor
python3 metax.py run-all --dry-run
python3 metax.py run-all
```

`run-all` chạy lần lượt và chạy tiếp được sau khi dừng. Thiếu công cụ hoặc database được ghi là BLOCKED. Trạng thái đó có nghĩa là bước chưa chạy, không phải bằng chứng rằng gene, taxon hoặc virus không có trong mẫu. Trên Mac, dùng `config/project_mac64.json`, tối đa 8 luồng và 40 GB cho mỗi tiến trình.

Đường dẫn để trống trong `reference_databases` sẽ tắt bước tương ứng.

## Mức phân tích

Mức nghiên cứu đủ cho thành phần quần xã, genome bin và chức năng mã hóa chính của một metagenome nước:

1. Kiểm tra reads, fastp, MEGAHIT, mapping, QUAST và Prodigal.
2. Phân loại trên reads bằng Kraken2 và Bracken, cùng SSU bằng phyloFlash.
3. Nếu bài báo dùng MAG: MetaBAT2, MaxBin2, CONCOCT, DAS Tool, CheckM2, ngưỡng 50% completeness và 10% contamination, rồi CoverM và dRep.
4. Một cơ sở ortholog, KOfam hoặc eggNOG, cùng dbCAN cho enzyme thủy phân carbohydrate.
5. Virus, khi đó là một câu hỏi của bài: geNomad và CheckV.
6. Tên loài trên bin chỉ khi bài báo ghi tên loài: GTDB-Tk. skani chỉ đưa ứng viên ANI, không phải tên loài GTDB.

Các bước còn lại trả lời thêm một câu hỏi hẹp hơn: phân loại thứ hai, sinh vật nhân thực, vật chủ của virus, cụm gene sinh tổng hợp, kháng kháng sinh, và độ phong phú pathway trên reads. Một nghiên cứu mẫu nước vẫn đầy đủ khi không chạy các bước đó.

DNA là tiềm năng được mã hóa. Kết quả không cho biết gene nào đang hoạt động vào ngày lấy mẫu.

## Phần cứng

Ubuntu, mức nghiên cứu: 64 GB RAM và ít nhất 1 TB trống. Ubuntu, mức đầy đủ: 64 GB RAM, trần 48 GB cho mỗi container Nextflow, và một tác vụ nặng tại một thời điểm. GTDB-Tk cần 128 GB RAM, hoặc 64 GB RAM cộng swap ít nhất 128 GB. Swap chậm hơn RAM. Nên chừa 4 TB trước khi cài đủ database của mức đầy đủ.

MacBook Pro 64 GB dùng profile `mac64`. Không chạy GTDB-Tk, chỉ mục Kraken 270 GB, iPHoP, hoặc DRAM2 với UniRef90 trên máy xách tay.

## Trích dẫn

Truong Nguyen, CUNY School of Medicine.
