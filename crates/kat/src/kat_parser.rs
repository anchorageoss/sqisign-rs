/// A single KAT test vector entry.
#[derive(Debug, Clone)]
pub struct KatEntry {
    pub count: usize,
    pub seed: Vec<u8>,
    pub mlen: usize,
    pub msg: Vec<u8>,
    pub pk: Vec<u8>,
    pub sk: Vec<u8>,
    pub smlen: usize,
    pub sm: Vec<u8>,
}

/// Parse a `.rsp` file into a list of KAT entries.
pub fn parse_rsp(contents: &str) -> Vec<KatEntry> {
    let mut entries = Vec::new();
    let mut current: Option<KatEntryBuilder> = None;

    for line in contents.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            if let Some(builder) = current.take() {
                if let Some(entry) = builder.build() {
                    entries.push(entry);
                }
            }
            continue;
        }

        let Some((key, value)) = line.split_once('=') else {
            continue;
        };
        let key = key.trim();
        let value = value.trim();

        if key == "count" {
            if let Some(builder) = current.take() {
                if let Some(entry) = builder.build() {
                    entries.push(entry);
                }
            }
            current = Some(KatEntryBuilder::new(value.parse().unwrap()));
            continue;
        }

        let Some(builder) = current.as_mut() else {
            continue;
        };

        match key {
            "seed" => builder.seed = Some(hex::decode(value).unwrap()),
            "mlen" => builder.mlen = Some(value.parse().unwrap()),
            "msg" => builder.msg = Some(hex::decode(value).unwrap()),
            "pk" => builder.pk = Some(hex::decode(value).unwrap()),
            "sk" => builder.sk = Some(hex::decode(value).unwrap()),
            "smlen" => builder.smlen = Some(value.parse().unwrap()),
            "sm" => builder.sm = Some(hex::decode(value).unwrap()),
            _ => {}
        }
    }

    if let Some(builder) = current.take() {
        if let Some(entry) = builder.build() {
            entries.push(entry);
        }
    }

    entries
}

#[derive(Default)]
struct KatEntryBuilder {
    count: usize,
    seed: Option<Vec<u8>>,
    mlen: Option<usize>,
    msg: Option<Vec<u8>>,
    pk: Option<Vec<u8>>,
    sk: Option<Vec<u8>>,
    smlen: Option<usize>,
    sm: Option<Vec<u8>>,
}

impl KatEntryBuilder {
    fn new(count: usize) -> Self {
        Self {
            count,
            ..Default::default()
        }
    }

    fn build(self) -> Option<KatEntry> {
        Some(KatEntry {
            count: self.count,
            seed: self.seed?,
            mlen: self.mlen?,
            msg: self.msg?,
            pk: self.pk?,
            sk: self.sk?,
            smlen: self.smlen?,
            sm: self.sm?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_rsp_entry() {
        let rsp = r#"# SQIsign_lvl1

count = 0
seed = 061550234D158C5EC95595FE04EF7A25767F2E24CC2BC479D09D86DC9ABCFDE7056A8C266F9EF97ED08541DBD2E1FFA1
mlen = 33
msg = D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8
pk = 07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B2029550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A010B
sk = 07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B2029550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A010B19943116DB5B4552B05B174969C61C9C8701000000000000000000000000000094F28A5533DF8872E3C7EFE3D45A175A0CFDFFFFFFFFFFFFFFFFFFFFFFFFFFFFF1959E3D67EADD79948DB766D9FFAF4D3FFDFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF0000000000000000000000000000000000000000000000000000000000000000358A8756E1CA2E31C2F3C879414AC08DF7EA0C1D732F9AE3D1AC4644E524340095A4F53D286FDE8A7226CE960C152C344888C963A457B02CAECA41C2672D76000365B548FB9C9E6C0E149BABA3EC7BC33B8F052B6B9D4F840A2AD67221C8F600464B9862D34ADF4D562F3836EBEFC4D8F874351B3E63A4DF9D33C0BBF9EB1800
smlen = 181
sm = 84228651F271B0F39F2F19F2E8718F31ED3365AC9E5CB303AFE663D0CFC11F0455D891B0CA6C7E653F9BA2667730BB77BEFE1B1A31828404284AF8FD7BAACC010001D974B5CA671FF65708D8B462A5A84A1443EE9B5FED7218767C9D85CEED04DB0A69A2F6EC3BE835B3B2624B9A0DF68837AD00BCACC27D1EC806A44840267471D86EFF3447018ADB0A6551EE8322AB30010202D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8
"#;
        let entries = parse_rsp(rsp);
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].count, 0);
        assert_eq!(entries[0].seed.len(), 48);
        assert_eq!(entries[0].mlen, 33);
        assert_eq!(entries[0].msg.len(), 33);
        assert_eq!(entries[0].pk.len(), 65);
        assert_eq!(entries[0].sk.len(), 353);
    }
}
