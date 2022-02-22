rule final_report:
    output: "cnv_sv/final_report.html"
    shell:
        "touch {output}"

