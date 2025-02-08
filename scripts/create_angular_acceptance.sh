#!/bin/bash

# スクリプトが存在するディレクトリを取得
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 相対パスを指定
DATA_DIR="${SCRIPT_DIR}/../../root"

"${SCRIPT_DIR}/../bin/acceptance_check" "${DATA_DIR}/for_acceptance_check_at735/mom735_flat_eta_lambda.root" 2212
"${SCRIPT_DIR}/../bin/acceptance_check" "${DATA_DIR}/for_acceptance_check_at735/mom735_flat_kp.root" 0
"${SCRIPT_DIR}/../bin/acceptance_check" "${DATA_DIR}/for_acceptance_check_at735/mom735_flat_k0n.root" 211
"${SCRIPT_DIR}/../bin/acceptance_check" "${DATA_DIR}/for_acceptance_check_at735/mom735_flat_pi+sigma-.root" 2112
"${SCRIPT_DIR}/../bin/acceptance_check" "${DATA_DIR}/for_acceptance_check_at735/mom735_flat_pi-sigma+.root" 2112
"${SCRIPT_DIR}/../bin/acceptance_check" "${DATA_DIR}/for_acceptance_check_at735/mom735_flat_pi-sigma+.root" 2212
"${SCRIPT_DIR}/../bin/acceptance_check" "${DATA_DIR}/for_acceptance_check_at735/mom735_flat_pi0sigma0.root" 2212
"${SCRIPT_DIR}/../bin/acceptance_check" "${DATA_DIR}/for_acceptance_check_at735/mom735_flat_pi0lambda.root" 2212