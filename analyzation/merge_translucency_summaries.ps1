$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$anaDir = Join-Path $root "analyzation"
$csvPath = Join-Path $anaDir "translucency_summary_merged.csv"
$xlsxPath = Join-Path $anaDir "translucency_summary_merged.xlsx"

New-Item -ItemType Directory -Path $anaDir -Force | Out-Null

# 1) Read all existing txt files under analyzation (including subfolders).
$txtFiles = Get-ChildItem -Path $anaDir -Recurse -File -Filter "*.txt" -ErrorAction SilentlyContinue |
    Where-Object { $_.FullName -ne $csvPath -and $_.FullName -ne $xlsxPath } |
    Sort-Object FullName

if (-not $txtFiles -or $txtFiles.Count -eq 0) {
    Write-Host "[WARN] No txt files found under: $anaDir"
    Write-Host "[HINT] Copy selected summary txt files into analyzation, then run analyzation_run.bat again."
    exit 0
}

# 2) Parse key-value records from each txt file.
$records = New-Object System.Collections.Generic.List[object]
$allKeys = New-Object System.Collections.Generic.List[string]
$keySet = @{}

function Add-Key([string]$k) {
    if ([string]::IsNullOrWhiteSpace($k)) { return }
    if (-not $keySet.ContainsKey($k)) {
        $keySet[$k] = $true
        $allKeys.Add($k) | Out-Null
    }
}

foreach ($file in $txtFiles) {
    $map = [ordered]@{}
    $map["FileName"] = $file.Name
    $map["FilePath"] = $file.FullName

    Add-Key "FileName"
    Add-Key "FilePath"

    $lines = Get-Content -Path $file.FullName -Encoding UTF8
    foreach ($raw in $lines) {
        $line = $raw.Trim()
        if ([string]::IsNullOrWhiteSpace($line)) { continue }

        # Pattern: [Final translucency edges] used X / Y edges.
        if ($line -match '^\[Final translucency edges\]\s*used\s*(\d+)\s*/\s*(\d+)\s*edges') {
            $map["FinalEdgesUsed"] = $Matches[1]
            $map["FinalEdgesTotal"] = $Matches[2]
            Add-Key "FinalEdgesUsed"
            Add-Key "FinalEdgesTotal"
            continue
        }

        # Generic key: value
        if ($line -match '^\s*([^:]+)\s*:\s*(.*)$') {
            $prefix = $Matches[1].Trim()
            $value = $Matches[2].Trim()

            # If value contains k=v pairs, flatten them.
            if ($value -match '=') {
                $parts = $value -split ','
                foreach ($p in $parts) {
                    $kv = $p.Trim()
                    if ($kv -match '^\s*([^=]+?)\s*=\s*(.*)\s*$') {
                        $k = ($prefix + "_" + $Matches[1].Trim()) -replace '\s+', '_'
                        $v = $Matches[2].Trim()
                        $map[$k] = $v
                        Add-Key $k
                    }
                }
            } else {
                $k = ($prefix) -replace '\s+', '_'
                $map[$k] = $value
                Add-Key $k
            }
            continue
        }
    }

    $records.Add([pscustomobject]$map) | Out-Null
}

if ($records.Count -eq 0) {
    Write-Host "[WARN] No parsable txt records in: $anaDir"
    exit 0
}

# 3) Export to CSV (Excel can open directly).
$records | Select-Object -Property $allKeys | Export-Csv -Path $csvPath -NoTypeInformation -Encoding UTF8
Write-Host "[OK] CSV written: $csvPath"

# 4) Try to create XLSX via Excel COM (if Excel exists).
try {
    $excel = New-Object -ComObject Excel.Application
    $excel.Visible = $false
    $excel.DisplayAlerts = $false
    $wb = $excel.Workbooks.Open($csvPath)
    $wb.SaveAs($xlsxPath, 51)  # 51 = xlOpenXMLWorkbook (*.xlsx)
    $wb.Close($true)
    $excel.Quit()
    [System.Runtime.Interopservices.Marshal]::ReleaseComObject($wb) | Out-Null
    [System.Runtime.Interopservices.Marshal]::ReleaseComObject($excel) | Out-Null
    Write-Host "[OK] XLSX written: $xlsxPath"
} catch {
    Write-Host "[WARN] Excel COM not available. Keep using CSV."
}

Write-Host "[DONE] Processed $($records.Count) txt files."
