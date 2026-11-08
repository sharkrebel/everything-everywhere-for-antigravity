---
name: deploy
description: Comprehensive step-by-step guide to deploy source code to GitHub and Ubuntu/Debian VPS with PostgreSQL, PM2, Nginx, SSL, UFW, and Fail2Ban.
---

# Quy trình Deploy Dự án - Toàn tập

Hướng dẫn này dành cho người mới, giúp bạn đưa dự án từ máy cá nhân lên VPS thông qua GitHub một cách an toàn và chuyên nghiệp.

---

## PHẦN 1: Quản lý Mã nguồn với GitHub

Việc sử dụng GitHub giúp bạn lưu trữ code an toàn và dễ dàng đồng bộ lên máy chủ VPS.

### 1.1 Tạo Repository trên GitHub (Chỉ làm lần đầu)
1. Truy cập [github.com/new](https://github.com/new).
2. Tên dự án: `tgv.bhxh`
3. Chế độ: **Private** (Riêng tư - Cực kỳ quan trọng để bảo mật code).
4. KHÔNG chọn bất kỳ tùy chọn khởi tạo nào (README, .gitignore, licenses).
5. Nhấn **Create repository**.

### 1.2 Kết nối Code local với GitHub (Chỉ làm lần đầu)
Mở terminal tại thư mục gốc của dự án (`c:\Users\ntruo\Desktop\Project\tgv.bhgiadinh`):

```powershell
# Khởi tạo git
git init

# Thêm địa chỉ server GitHub (Thay USERNAME bằng tên GitHub của bạn)
git remote add origin https://github.com/USERNAME/tgv.bhxh.git

# Đặt tên nhánh chính là main
git branch -M main

# Lưu toàn bộ code và đẩy lên GitHub
git add .
git commit -m "Initial commit: Khởi tạo dự án"
git push -u origin main
```

### 1.3 Cập nhật code sau này
Mỗi khi bạn sửa code xong và muốn đưa lên máy chủ:

```powershell
git add .
git commit -m "Mô tả nội dung bạn vừa sửa (ví dụ: Sửa lỗi hiển thị)"
git push
```

---

## PHẦN 2: Thiết lập Máy chủ VPS (Ubuntu/Debian)

### 2.1 Cài đặt môi trường cần thiết
Đăng nhập vào VPS qua SSH và chạy các lệnh sau:

```bash
# Cập nhật hệ thống
sudo apt update && sudo apt upgrade -y

# Cài đặt công cụ hỗ trợ biên dịch (Cần cho bcrypt và các thư viện hệ thống)
sudo apt install -y build-essential python3

# Cài đặt Node.js 20 (Dùng để chạy app)
curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
sudo apt install -y nodejs

# Cài đặt PostgreSQL (Cơ sở dữ liệu)
sudo apt install -y postgresql postgresql-contrib

# Cài đặt Nginx (Giúp truy cập web qua domain)
sudo apt install -y nginx

# Cài đặt LibreOffice (Dùng để chuyển đổi Docx sang PDF)
sudo apt install -y libreoffice-writer libreoffice-java-common

# Cài đặt Git (Để tải code từ GitHub)
sudo apt install -y git

# Cài đặt PM2 (Giúp app tự động chạy lại nếu server bị sập)
sudo npm install -g pm2
```

### 2.2 Cấu hình Cơ sở dữ liệu (PostgreSQL)
```bash
sudo -u postgres psql

# Chạy các lệnh SQL sau (Thay đổi 'mat_khau_cua_ban' thành mật khẩu bảo mật)
CREATE USER bhgiadinh WITH PASSWORD 'mat_khau_cua_ban';
CREATE DATABASE bhgiadinh_db OWNER bhgiadinh;
\q
```

---

## PHẦN 3: Triển khai Ứng dụng lần đầu

### 3.1 Tải Code về VPS
```bash
cd /var/www
# Thay USERNAME bằng tên GitHub của bạn
git clone https://github.com/USERNAME/tgv.bhxh.git
cd tgv.bhxh
```

### 3.2 Cấu hình Biến môi trường (.env) cho Backend
Tại thư mục gốc dự án:
```bash
cp .env.example .env
nano .env
```
**Chỉnh sửa nội dung cho đúng với VPS:**
```env
DATABASE_URL="postgresql://bhgiadinh:mat_khau_cua_ban@localhost:5432/bhgiadinh_db"
JWT_SECRET="tạo_mã_bí_mật_của_bạn_ít_nhất_32_ký_tự"

# Tạo khóa mã hóa file (Dùng lệnh bên dưới để tạo 64 ký tự hex)
FILE_ENCRYPTION_KEY="64_ký_tự_hex_vừa_tạo"

PORT=3001
FRONTEND_URL="https://ten-mien-cua-ban.vn"
NODE_ENV=production
```

# Bảo mật file .env (Chỉ cho phép chủ sở hữu đọc/ghi)
chmod 600 .env
> **Mẹo tạo khóa mã hóa:** Chạy lệnh `node -e "console.log(require('crypto').randomBytes(32).toString('hex'))"` để lấy mã.

### 3.3 Build Backend & Database
```bash
npm install --legacy-peer-deps
npx prisma generate
npx prisma migrate deploy

# Chạy seed dữ liệu (Chỉ làm lần đầu)
npx ts-node prisma/seed.ts

# Cập nhật danh sách bệnh viện
node prisma/update-hospitals.js

# Biên dịch code
npm run build

# --- TẠO THƯ MỤC STORAGE & CẤP QUYỀN (QUAN TRỌNG) ---
mkdir -p /var/www/tgv.bhxh/storage
sudo chown -R $USER:$USER /var/www/tgv.bhxh/storage
chmod -R 755 /var/www/tgv.bhxh/storage

# Chạy với PM2 (Chỉ nghe từ localhost để bảo mật)
pm2 start dist/src/main.js --name "bhgiadinh-api"
```

### 3.4 Cấu hình & Build Frontend
```bash
cd frontend
cp .env.production.example .env.production
nano .env.production
```
**Nội dung `.env.production`:**
```env
NEXT_PUBLIC_API_URL=https://ten-mien-cua-ban.vn/api
BACKEND_URL=http://localhost:3001
NEXT_PUBLIC_SITE_URL=https://ten-mien-cua-ban.vn
```

# Bảo mật file .env (Chỉ cho phép chủ sở hữu đọc/ghi)
chmod 600 .env.production

**Cài đặt và chạy:**
```bash
npm install --legacy-peer-deps
npm run build

# Chạy với PM2 (Binding vào localhost để bảo mật)
pm2 start "npm run start -- -H 127.0.0.1" --name "bhgiadinh-web"
```

---

## PHẦN 4: Cấu hình Tên miền (Nginx)

```bash
sudo nano /etc/nginx/sites-available/bhgiadinh
```
**Dán nội dung sau vào (Thay `ten-mien-cua-ban.vn` bằng tên miền thật):**
> [!IMPORTANT]
> **Lưu ý cực kỳ quan trọng:** Khi copy nội dung bên dưới, bạn **CHỈ COPY** phần code từ dòng `server {` đến hết dấu `}`. **KHÔNG COPY** các dấu gạch chéo ngược ( \` \` \` nginx ) ở đầu và cuối, nếu không Nginx sẽ báo lỗi không chạy được.

```nginx
# Đặt trong http block (nginx.conf)
limit_req_zone $binary_remote_addr zone=registration_limit:10m rate=10r/m;

# Đặt trọng nginx/sites-available/bhgiadinh
server {
    server_name muabhyt.vn;
    server_tokens off;
    
    # Security headers
    add_header X-Frame-Options "SAMEORIGIN" always;
    add_header X-Content-Type-Options "nosniff" always;
    add_header X-XSS-Protection "1; mode=block" always;
    
    # Timeouts
    client_body_timeout 12;
    client_header_timeout 12;
    
    # File upload
    client_max_body_size 5M;
    
    # Hide .env
    location ~ /\.env {
        deny all;
    }
    
    # Rate limited BHYT/BHXH endpoints - ĐẶT TRƯỚC /api
    location ~ ^/api/orders/(bhyt|bhxh) {
        limit_req zone=registration_limit burst=15 nodelay;
        limit_req_status 429;
        
        proxy_pass http://localhost:3001;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_cache_bypass $http_upgrade;
    }
    
    # Backend API
    location /api {
        proxy_pass http://localhost:3001;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }
    
    # Frontend (Next.js)
    location / {
        proxy_pass http://localhost:3000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_cache_bypass $http_upgrade;
    }

    listen 443 ssl http2;
    ssl_certificate /etc/letsencrypt/live/muabhyt.vn/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/muabhyt.vn/privkey.pem;
    include /etc/letsencrypt/options-ssl-nginx.conf;
    ssl_dhparam /etc/letsencrypt/ssl-dhparams.pem;
}

server {
    listen 80;
    server_name muabhyt.vn;
    return 301 https://$host$request_uri;
}
```

```bash
# Kích hoạt cấu hình lần đầu - 3 lệnh
# Update file chỉ cần nhập 2 lệnh cuối
sudo ln -s /etc/nginx/sites-available/bhgiadinh /etc/nginx/sites-enabled/ 
sudo nginx -t
sudo systemctl reload nginx

```

---

## PHẦN 5: Cấu hình PM2 tự khởi động lại

Để đảm bảo các ứng dụng chạy bằng PM2 tự động khởi động lại sau khi VPS reboot:

```bash
# Lưu trạng thái PM2 (Để khi reboot server app tự chạy lại)
pm2 save
pm2 startup
# (Lưu ý: Sau lệnh startup, PM2 sẽ hiện ra 1 dòng lệnh, hãy copy và chạy dòng lệnh đó)
```

---

## PHẦN 6: Cấu hình SSL (HTTPS) - Khuyên dùng

Để website chuyên nghiệp và bảo mật (có biểu tượng ổ khóa), bạn nên cài đặt SSL miễn phí bằng Certbot.

```bash
# Cài đặt Certbot cho Nginx
sudo apt install certbot python3-certbot-nginx -y

# Chạy lệnh cấp chứng chỉ (Thay ten-mien-cua-ban.vn bằng tên miền thật)
# Lệnh này sẽ tự động sửa file cấu hình Nginx cho bạn
sudo certbot --nginx -d ten-mien-cua-ban.vn
```
> **Lưu ý:** Khi chạy lệnh trên, Certbot sẽ hỏi email (để báo tin sắp hết hạn) và hỏi bạn có muốn tự chuyển hướng HTTP sang HTTPS không (nên chọn **2 - Redirect**).

---

## PHẦN 7: Bảo mật với Tường lửa (UFW) - Khuyên dùng

Để bảo vệ máy chủ khỏi các truy cập trái phép, bạn nên bật tường lửa và chỉ mở các cổng cần thiết.

```bash
# Cài đặt UFW (nếu chưa có)
sudo apt update
sudo apt install ufw -y

# Thiết lập mặc định
sudo ufw default deny incoming
sudo ufw default allow outgoing

# Mở cổng SSH (Cực kỳ quan trọng để không bị khóa máy)
sudo ufw allow 24700/tcp

# Mở cổng cho Web
sudo ufw allow 80/tcp
sudo ufw allow 443/tcp

# Kích hoạt tường lửa (Nhấn 'y' khi được hỏi)
sudo ufw enable

# Kiểm tra trạng thái
sudo ufw status
```

---

## PHẦN 8: Setup Fail2Ban (Tự động block IP tấn công)

Fail2Ban giúp tự động hóa việc chặn các địa chỉ IP có dấu hiệu tấn công (như dò mật khẩu SSH hoặc spam request).

```bash
# 8.1 Cài đặt
sudo apt install fail2ban -y

# 8.2 Tạo file cấu hình riêng
sudo nano /etc/fail2ban/jail.local
```

**Dán nội dung cấu hình sau vào:**

```ini
[DEFAULT]
bantime = 3600       # Block IP trong 1 giờ
findtime = 600       # Trong vòng 10 phút
maxretry = 5         # Nếu vi phạm 5 lần thì block

[sshd]
enabled = true
port = 24700         # QUAN TRỌNG: Khớp với cổng SSH đã đổi ở Phần 7
logpath = /var/log/auth.log

[nginx-http-auth]
enabled = true
filter = nginx-http-auth
port = http,https
logpath = /var/log/nginx/error.log

[nginx-limit-req]
enabled = true
filter = nginx-limit-req
port = http,https
logpath = /var/log/nginx/error.log
maxretry = 3

[nginx-botsearch]
enabled = true
port = http,https
filter = nginx-botsearch
logpath = /var/log/nginx/access.log
maxretry = 2
```

```bash
# 8.3 Khởi động và kích hoạt
sudo systemctl enable fail2ban
sudo systemctl start fail2ban

# 8.4 Kiểm tra trạng thái
sudo fail2ban-client status
sudo fail2ban-client status sshd
```

---

## PHẦN 9: Cập nhật dự án (Về sau)

Mỗi khi có code mới trên GitHub, hãy chạy các lệnh này trên VPS:

```bash
cd /var/www/tgv.bhxh
git pull origin main
git reset --hard (nếu pull conflict)

# --- QUAN TRỌNG: Nếu bạn có thêm thư viện mới (như throttler, v.v.) Bạn BẮT BUỘC phải chạy lệnh install và build lại:
npm install --legacy-peer-deps
# Cập nhật Backend
npx prisma generate
npx prisma migrate deploy
npx prisma db push (nếu chưa có migrate)
npm run build 
pm2 restart bhgiadinh-api

# Cập nhật Frontend
cd frontend
npm install --legacy-peer-deps
npm run build
pm2 restart bhgiadinh-web

# Warmup Sitemap Shards (Tránh cache lạnh sau deploy)
curl -s http://localhost:3000/sitemap.xml > /dev/null &
curl -s http://localhost:3000/sitemap/pages.xml > /dev/null &
curl -s http://localhost:3000/sitemap/provinces.xml > /dev/null &
curl -s http://localhost:3000/sitemap/districts.xml > /dev/null &
curl -s http://localhost:3000/sitemap/wards-1.xml > /dev/null &
curl -s http://localhost:3000/sitemap/wards-2.xml > /dev/null &
curl -s http://localhost:3000/sitemap/wards-3.xml > /dev/null &

pm2 save

# --- TUỲ CHỌN: Cập nhật "Sạch" (Clean Deploy) ---
# Nếu web không cập nhật đúng code mới, hãy dùng cách này:
cd /var/www/tgv.bhxh/frontend
rm -rf .next
npm run build
pm2 restart bhgiadinh-web
# Warmup sitemap shards sau clean deploy:
curl -s http://localhost:3000/sitemap.xml > /dev/null &
curl -s http://localhost:3000/sitemap/pages.xml > /dev/null &
curl -s http://localhost:3000/sitemap/provinces.xml > /dev/null &
curl -s http://localhost:3000/sitemap/districts.xml > /dev/null &
curl -s http://localhost:3000/sitemap/wards-1.xml > /dev/null &
curl -s http://localhost:3000/sitemap/wards-2.xml > /dev/null &
curl -s http://localhost:3000/sitemap/wards-3.xml > /dev/null &

# --- TUỲ CHỌN: Quản lý tập trung với Ecosystem ---
# Dùng lệnh này để khởi động/restart cả 2 app cùng lúc:
cd /var/www/tgv.bhxh
pm2 start ecosystem.config.js --env production
# Restart lại cả 2:
pm2 restart ecosystem.config.js --env production
pm2 save

```
> [!TIP]
> **Quy trình chuẩn khi thêm tính năng mới:**
> 1. Tại máy cá nhân: Cài đặt thư viện -> Viết code -> `git push`.
> 2. Tại VPS: `git pull` -> `npm install` -> `npm run build` -> `pm2 restart`.

---

### CÁC LƯU Ý QUAN TRỌNG & LỖI THƯỜNG GẶP

1. **Cảnh báo `sudo: unable to resolve host...`:**
   - **Hiện tượng:** Khi chạy lệnh `sudo`, bạn thấy dòng chữ `unable to resolve host [tên_vps]`.
   - **Xử lý:** Đây chỉ là cảnh báo về cấu hình tên máy (hostname) của VPS, **không ảnh hưởng** đến việc cài đặt. Bạn có thể bỏ qua hoặc sửa bằng cách chạy lệnh: `echo "127.0.0.1 $(hostname)" | sudo tee -a /etc/hosts`.

2. **Bảo mật .env:** Tuyệt đối không bao giờ push tệp `.env` lên GitHub.
3. **Backup DB:** Thường xuyên sao lưu PostgreSQL bằng lệnh `pg_dump`.
4. **Log lỗi:** Nếu app không chạy, xem log bằng:
   - `pm2 logs bhgiadinh-api`
   - `pm2 logs bhgiadinh-web`
5. **Bộ nhớ RAM:** Dự án cần ít nhất 2GB RAM để quá trình `npm run build` không bị treo.
6. **Thông báo cập nhật `npm`:**
   - Khi chạy lệnh, bạn có thể thấy `npm notice New major version of npm available!`.
   - **Xử lý:** Bạn **không cần cập nhật**. Phiên bản đi kèm với Node.js 20 đã rất ổn định. Việc cập nhật lên bản mới nhất (ví dụ bản 11) đôi khi gây ra một số thay đổi không cần thiết trên máy chủ production.

8. **Lỗi kết nối PostgreSQL (`psql: error...`):**
   - **Hiện tượng:** Khi chạy `sudo -u postgres psql`, bạn thấy báo lỗi `No such file or directory` liên quan đến `5432`.
   - **Xử lý:** Có thể dịch vụ Postgres chưa khởi động. Hãy chạy:
     ```bash
     sudo systemctl start postgresql
     sudo systemctl enable postgresql
     ```
   - **Lưu ý:** Nếu thấy cảnh báo `could not change directory to "/root": Permission denied`, bạn có thể **bỏ qua**. Đó chỉ là thông báo Postgres không có quyền vào thư mục `root` của bạn, nhưng nó vẫn sẽ vào được trình quản lý DB.

10. **Lỗi `npm error ERESOLVE`:**
    - **Hiện tượng:** Khi chạy `npm install`, hệ thống báo lỗi xung đột phiên bản (thường là giữa `nodemailer` và `mailtrap`).
    - **Xử lý:** Đây là do bản npm mới kiểm tra phụ thuộc rất khắt khe. Bạn chỉ cần thêm tham số `--legacy-peer-deps` vào lệnh cài đặt:
      ```bash
      npm install --legacy-peer-deps
      ```

13. **Cấu hình HTTPS:** Luôn luôn nên thực hiện bước ở **PHẦN 6** sau khi cấu hình xong Nginx để đảm bảo an toàn dữ liệu khách hàng.

14. **Về các thư viện như `pdfkit`, `exceljs`, `pdf-lib`:**
    - Bạn **KHÔNG CẦN** cài đặt chúng thủ công bằng `apt`.
    - Đây là các thư viện Node.js, chúng sẽ tự động được cài đặt vào thư mục `node_modules` khi bạn chạy lệnh `npm install --legacy-peer-deps`.
    - Chỉ những phần mềm hệ thống (như Postgres, Nginx, LibreOffice) mới cần cài bằng `sudo apt`.

15. **Lỗi `ERR_SSL_VERSION_OR_CIPHER_MISMATCH` khi dùng Cloudflare:**
    - **Hiện tượng:** Truy cập web qua Cloudflare báo lỗi về SSL/Cipher.
    - **Nguyên nhân:** Cloudflare đang trong quá trình khởi tạo Chứng chỉ SSL (thường mất 5-15 phút sau khi đổi Nameserver).
    - **Xử lý:** 
        - Đợi khoảng 15 phút để Cloudflare hoàn tất việc cấp phát "Universal SSL".
        - Trong bảng điều khiển Cloudflare, vào mục **SSL/TLS** -> **Overview**, chọn chế độ **Full (Strict)** nếu bạn đã chạy Certbot trên VPS, hoặc **Flexible** nếu chưa chạy Certbot. Khuyên dùng **Full (Strict)** để bảo mật cao nhất.

16. **Nâng cấp Rate Limiting (Chống Spam/Brute Force):**
    - Nếu sau này bạn thấy web bị spam hoặc có nhiều yêu cầu lạ, hãy cài đặt thêm `@nestjs/throttler`.
    - **Quy trình:** 
        1. Máy cá nhân: `npm install @nestjs/throttler`.
        2. Config trong `app.module.ts`: Khai báo `ThrottlerModule.forRoot({...})`.
        3. Kích hoạt trong `main.ts` hoặc từng Controller: Dùng `ThrottlerGuard`.
        4. Deploy: Chạy theo đúng **PHẦN 9** (git pull -> install -> build -> restart).

17. **Cấu hình Bảo mật Nâng cao (HSTS & CSP):**
    - **HSTS:** Nên bật trên **Cloudflare** (SSL/TLS -> Edge Certificates -> Always Use HTTPS & HSTS) thay vì cấu hình trong code để đạt hiệu quả tốt nhất.
    - **CSP:** Giúp chống XSS. Nếu muốn thêm, cần liệt kê đủ các nguồn (nhà cung cấp font, icon, biểu đồ) để tránh làm hỏng giao diện.

18. **Bảo mật Server mức độ cao (SSH Key & Fail2Ban):**
    - **SSH Key:** Thay vì dùng mật khẩu, hãy dùng chìa khóa số (Private/Public Key). Sau khi cài xong, hãy cấu hình `PasswordAuthentication no` để chặn hoàn toàn mọi nỗ lực dò mật khẩu.
    - **Fail2Ban:** Công cụ này giống như "Bảo vệ tự động". Nó sẽ theo dõi log hệ thống, nếu thấy IP nào nhập sai pass 3-5 lần, nó sẽ tự động khóa (ban) IP đó trong một khoảng thời gian. Kết hợp với việc bạn đã đổi cổng SSH (24700), server sẽ cực kỳ an toàn.

19. **Khắc phục lỗ hổng từ Pen-test:**
    - **Lộ thông tin phiên bản:** Đã xóa `name` và `version` tại API gốc (`/api`).
    - **X-Powered-By & Security Headers (Backend):** Đã cài đặt và kích hoạt **Helmet** trong `main.ts`.
    - **Security Headers (Frontend):** Đã bổ sung bộ 3 tiêu đề bảo mật cơ bản (`X-Frame-Options`, v.v.) vào `next.config.js`.

---

## PHẦN 10: Bảo trì và Cập nhật định kỳ

Dự án cần được bảo trì để đảm bảo an toàn và hiệu suất.

### 10.1 Kiểm tra lỗi bảo mật (NPM Audit)
Thỉnh thoảng bạn nên kiểm tra xem các thư viện mình dùng có bị lỗ hổng bảo mật nào mới phát hiện không:
```bash
# Chỉ kiểm tra (Nên làm thường xuyên)
npm audit

# Tự động sửa lỗi (Cẩn thận: có thể làm lỗi code)
# Chỉ nên chạy ở máy cá nhân, test kỹ trước khi push lên VPS
npm audit fix
```

### 10.2 Cập nhật Node.js
Node.js phiên bản 20 là bản ổn định lâu dài (LTS). Bạn không cần cập nhật hàng tuần, nhưng khoảng 3-6 tháng nên chạy lệnh sau trên VPS một lần:
```bash
sudo apt update
sudo apt upgrade nodejs -y
```

### 10.3 Cập nhật Next.js & Thư viện
Việc cập nhật các bản lớn (Major version) của Next.js (ví dụ từ 14 lên 15) thường có lỗi vỡ giao diện.
*   **Quy trình:** Chỉ cập nhật tại **máy cá nhân** trước (`npm install next@latest`), sửa hết lỗi đỏ, sau đó mới `git push` lên VPS.
*   **KHÔNG** chạy cập nhật phiên bản lớn trực tiếp trên VPS.

---

## PHẦN 11: Giám sát Hệ thống (Monitoring)

Để đảm bảo website luôn hoạt động ổn định và phát hiện sớm các truy cập bất thường (tấn công), bạn nên thiết lập các công cụ giám sát sau:

### 11.1 Quản lý Log tự động (PM2 Logrotate)
PM2 ghi lại toàn bộ nhật ký ứng dụng. Nếu không kiểm soát, các file log này sẽ ngày càng lớn và làm đầy ổ cứng VPS.
```bash
# Cài đặt module tự động nén và xoay vòng log
pm2 install pm2-logrotate

# Cấu hình: Giữ tối đa 10 file log, mỗi file tối đa 10MB
pm2 set pm2-logrotate:max_size 10M
pm2 set pm2-logrotate:retain 10
```

### 11.2 Kiểm tra nhanh truy cập bất thường (Log Analysis)

Để theo dõi và xử lý các hành vi nghi vấn, bạn có thể dùng các lệnh sau:

```bash
# 1. Xem nhật ký truy cập theo thời gian thực (Real-time)
sudo tail -f /var/log/nginx/access.log

# 2. Xem nhật ký lỗi (Error log) - Cực kỳ hữu ích khi web báo lỗi 502/504
sudo tail -f /var/log/nginx/error.log

# 3. Lọc danh sách 10 IP truy cập nhiều nhất (Dễ nhận diện Bot/Spam)
sudo awk '{print $1}' /var/log/nginx/access.log | sort | uniq -c | sort -nr | head -10

# 4. Kiểm tra mọi hoạt động của một IP cụ thể (Thay IP_CẦN_KIỂM_TRA)
sudo grep "IP_CẦN_KIỂM_TRA" /var/log/nginx/access.log

# 5. Tạo lệnh tắt (Alias) để soi nhanh các lỗi bảo mật
echo "alias check-suspicious='sudo grep -E \"(404|403|429)\" /var/log/nginx/access.log | tail -50'" >> ~/.bashrc
source ~/.bashrc
```

### 11.3 Dashboard xem trực quan (GoAccess)
GoAccess giúp bạn xem biểu đồ thống kê truy cập ngay trên giao diện web hoặc terminal rất đẹp mắt.
```bash
# Cài đặt
sudo apt install goaccess -y

# Xem báo cáo trực tiếp từ Nginx log
sudo goaccess /var/log/nginx/access.log --log-format=COMBINED
```

### 11.4 Nhận thông báo qua Telegram (Khuyên dùng về sau)
Bạn có thể cài đặt để server tự nhắn tin vào điện thoại mỗi khi có lỗi nghiêm trọng (ví dụ: Database bị sập). 
*   **Cách làm:** Viết một script nhỏ gọi API của Telegram Bot và đặt vào Cronjob hoặc tích hợp vào code Backend.

---

## PHẦN 12: Bảo mật SSH nâng cao (SSH Key)

Đây là bước bảo mật mạnh nhất cho máy chủ, giúp thay thế mật khẩu bằng chìa khóa số vô cùng bảo mật.

### 12.1 Tạo và Copy SSH Key (Làm tại Máy cá nhân)
Nếu bạn chưa có SSH Key, hãy mở Terminal/PowerShell ở máy mình:
```powershell
# Tạo cặp khóa (Nhấn Enter liên tục)
ssh-keygen -t ed25519

# Copy khóa lên VPS (Thay USER và IP_VPS)
# Nếu dùng Windows:
type $env:USERPROFILE\.ssh\id_ed25519.pub | ssh -p 24700 USER@IP_VPS "mkdir -p ~/.ssh && cat >> ~/.ssh/authorized_keys"
```

### 12.2 Cấu hình Hardening (Làm trên VPS)
Sau khi đã chắc chắn đăng nhập được bằng Key (không cần gõ pass), hãy khóa server lại:
```bash
sudo nano /etc/ssh/sshd_config
```

**Tìm và sửa/thêm các dòng sau:**
```ini
# 1. Tắt đăng nhập bằng mật khẩu
PasswordAuthentication no

# 2. Bật đăng nhập bằng Public Key
PubkeyAuthentication yes

# 3. Tắt các phương thức xác thực kém an toàn
ChallengeResponseAuthentication no
UsePAM no

# 4. Tắt Root Login (Khuyên dùng cực cao - Đảm bảo bạn có user sudo)
PermitRootLogin no

# 5. Chỉ cho phép các user nhất định (Optional)
# AllowUsers ten_user_cua_ban
```

**Áp dụng cấu hình:**
> [!CAUTION]
> **CẢNH BÁO:** Đừng tắt Session hiện tại! Hãy mở một Terminal mới và thử `ssh` vào lại. Nếu vào được mới được đóng Terminal cũ.

```bash
sudo systemctl restart ssh
```
