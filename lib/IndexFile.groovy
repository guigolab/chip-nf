class IndexFile {

  static parse(def f, def filter) {
    def m = [:]
    f.readLines().each {
        def (path, meta) = it.split('\t')
        def metaDict = [:]
        meta.split('; ').each {
            def (k, v) = it.split('=')
            metaDict[k] = v
        }
        if ( ! filterLine(metaDict, filter) ) {
          m[path] = metaDict
        }
    }
    return m.entrySet()
  }

  static format(def files) {
    def out = files.collect { path, meta ->
        formatLine(path, meta)
    }
    return out.join('\n')
  }

  static formatLine(def path, def meta) {
    return "${path}\t${meta.entrySet().join('; ')};"
  }

  static filterLine(def meta, def filter) {
    def res = filter.collect { k,v ->
      if (k in meta) {
        if ( v instanceof String )
          meta[k] ==~ /${v}/
        else if ( v instanceof List )
          true in v.collect { meta[k] ==~ /${it}/ }
      }
    }
    return false in res
  }

}
