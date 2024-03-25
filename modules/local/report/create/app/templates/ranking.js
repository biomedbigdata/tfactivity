const initRanking = function (primary_ranking, secondary_ranking) {
    const assayChips = Array.from(document.querySelectorAll('[id^="assay-"]'));

    const updatePrimaryRanking = async function (activeAssays) {
        const primaryScores = Object.entries(primary_ranking).reduce(function (acc, [primary, dcgs]) {
            const score = activeAssays.reduce(function (acc, assay) {
                return acc + (dcgs[assay] || 0);
            }, 0) / activeAssays.length;

            acc[primary] = score;
            return acc;
        }, {});

        const primaryOrder = Object.keys(primaryScores).sort(function (a, b) {
            return primaryScores[b] - primaryScores[a];
        });

        const primaryRank = primaryOrder.reduce(function (acc, primary, index) {
            acc[primary] = index;
            return acc;
        }, {});

        Object.entries(primaryRank).forEach(function ([primary, rank]) {
            const primaryCard = document.getElementById(`primary-${primary}`);
            primaryCard.style.order = rank;
        });
    };

    const updateSecondaryRanking = async function (activeAssays) {
        for (const primary in secondary_ranking) {
            const currentRanking = secondary_ranking[primary];
            const secondaryScores = Object.entries(currentRanking).reduce(function (acc, [secondary, dcgs]) {
                const score = activeAssays.reduce(function (acc, assay) {
                    return acc + (dcgs[assay] || 0);
                }, 0) / activeAssays.length;

                acc[secondary] = score;
                return acc;
            }, {});

            const tgOrder = Object.keys(secondaryScores).sort(function (a, b) {
                return secondaryScores[b] - secondaryScores[a];
            });

            const tgRank = tgOrder.reduce(function (acc, secondary, index) {
                acc[secondary] = index;
                return acc;
            }, {});

            const showedSecondaries = []

            Object.entries(tgRank).forEach(function ([secondary, rank]) {
                const tgChip = document.getElementById(`primary-${primary}-secondary-${secondary}`);
                tgChip.style.order = rank;
                tgChip.style.display = rank <= 30 ? 'inline-block' : 'none';
                if (rank <= 30) {
                    showedSecondaries.push(secondary)
                }
            });

            const gProfilerLink = `https://biit.cs.ut.ee/gprofiler/gost?query=${showedSecondaries.join('%0D')}`;
            const gProfilerLinkElement = document.getElementById(`primary-${primary}-gprofiler`);
            gProfilerLinkElement.href = gProfilerLink;
        }
    }

    const allAssays = assayChips.map(function (chip) {
        return chip.textContent;
    });

    updatePrimaryRanking(allAssays);
    updateSecondaryRanking(allAssays);

    assayChips.forEach(function (element) {
        element.addEventListener('click', function () {
            if (element.classList.contains('c-not-allowed')) {
                return;
            }

            element.classList.toggle('active');

            const activeAssayChips = assayChips.filter(function (chip) {
                return chip.classList.contains('active');
            });

            if (activeAssayChips.length === 1) {
                activeAssayChips[0].classList.remove('c-hand');
                activeAssayChips[0].classList.add('c-not-allowed');
            } else {
                activeAssayChips.forEach(function (chip) {
                    chip.classList.remove('c-not-allowed');
                    chip.classList.add('c-hand');
                });
            }

            const activeAssays = activeAssayChips.map(function (chip) {
                return chip.textContent;
            });

            updatePrimaryRanking(activeAssays);
            updateSecondaryRanking(activeAssays);
        });
    });

    const filterPrimary = document.getElementById('filter-primary');
    filterPrimary.addEventListener('input', function () {
        const filter = filterPrimary.value.toLowerCase();
        const primaryCards = Array.from(document.querySelectorAll('.primary-card'));

        primaryCards.forEach(function (card) {
            const primary = card.id.replace('primary-', '');
            if (primary.toLowerCase().includes(filter)) {
                card.style.display = 'block';
            } else {
                card.style.display = 'none';
            }
        });
    });
};