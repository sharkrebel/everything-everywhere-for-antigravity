# Script cài đặt / đồng bộ skills sang Claude Code (~/.claude/skills/)
param (
    [string]$TargetDir = "$env:USERPROFILE\.claude\skills",
    [string]$ClaudeConfig = "$env:USERPROFILE\.claude"
)

Write-Host "==========================================================" -ForegroundColor Cyan
Write-Host " [Claude Code] Đồng bộ bộ kỹ năng AI Coding Skills Kit" -ForegroundColor Green
Write-Host "==========================================================" -ForegroundColor Cyan

if (-not (Test-Path $TargetDir)) {
    New-Item -ItemType Directory -Path $TargetDir -Force | Out-Null
}

$CurrentDir = Split-Path -Parent $PSScriptRoot
$SkillsSource = Join-Path $CurrentDir "skills"
$RulesSource = Join-Path $CurrentDir "rules"
$AgentsSource = Join-Path $CurrentDir "agents"

Write-Host ">> Nguồn skills: $SkillsSource"
Write-Host ">> Đích Claude: $TargetDir"

$skills = Get-ChildItem -Path $SkillsSource -Directory
$count = 0
foreach ($s in $skills) {
    $dest = Join-Path $TargetDir $s.Name
    Copy-Item -Path $s.FullName -Destination $dest -Recurse -Force
    $count++
}

Write-Host ">> Đã đồng bộ thành công $count skills sang Claude Code!" -ForegroundColor Green

# Đồng bộ AGENTS.md và SOUL.md vào ~/.claude/
if (Test-Path "$RulesSource\AGENTS.md") {
    Copy-Item -Path "$RulesSource\AGENTS.md" -Destination "$ClaudeConfig\AGENTS.md" -Force
    Copy-Item -Path "$RulesSource\AGENTS.md" -Destination "$ClaudeConfig\CLAUDE.md" -Force
    Write-Host ">> Đã cập nhật rules AGENTS.md & CLAUDE.md!" -ForegroundColor Yellow
}

if (Test-Path "$RulesSource\SOUL.md") {
    Copy-Item -Path "$RulesSource\SOUL.md" -Destination "$ClaudeConfig\SOUL.md" -Force
    Write-Host ">> Đã cập nhật SOUL.md!" -ForegroundColor Yellow
}

# Đồng bộ agents vào ~/.claude/agents/
if (Test-Path $AgentsSource) {
    $targetAgents = "$ClaudeConfig\agents"
    if (-not (Test-Path $targetAgents)) {
        New-Item -ItemType Directory -Path $targetAgents -Force | Out-Null
    }
    $agents = Get-ChildItem -Path $AgentsSource -Filter "*.md"
    foreach ($a in $agents) {
        Copy-Item -Path $a.FullName -Destination "$targetAgents\$($a.Name)" -Force
    }
    Write-Host ">> Đã cập nhật $($agents.Count) agent personas sang Claude Code!" -ForegroundColor Yellow
}

Write-Host " Hoàn tất cấu hình Claude Code!" -ForegroundColor Cyan
