static void half_convert_reference(unsigned logn, hw *d, const uint32_t *f, unsigned len, unsigned scale)
{
    unsigned n = 1U << logn, sch = (scale + 13) / 31, bits = (scale + 13) % 31;
    for (unsigned i = 0; i < n; i++) {
        uint32_t w[6] = {0}, sign = -(f[(len - 1) * n + i] >> 30) >> 1;
        for (unsigned j = 0; j < len; j++) {
            uint32_t x = f[j * n + i];
            for (unsigned k = 0; k < 6; k++)
                w[k] |= x & (uint32_t)hw_mask(j + 3 == sch + k);
        }
        for (unsigned k = 0; k < 6; k++)
            w[k] |= sign & (uint32_t)hw_mask(sch + k >= len + 3);
        hw lo = (hw)w[0] | (hw)w[1] << 31 | (hw)w[2] << 62 | (hw)w[3] << 93 | (hw)w[4] << 124;
        uint64_t hi = ((uint64_t)w[4] >> 4) | (uint64_t)w[5] << 27;
        d[i] = (lo >> bits) | (((hw)hi << 1) << (127 - bits));
    }
}
